import ../ttf, os, nimBMP, streams, unicode, sequtils

proc findFileInPaths(name: string, paths: varargs[string]): string =
    for p in paths:
        result = p / name
        if fileExists(result): return
    result = ""

proc findFile(name: string): string = findFileInPaths(name, ".", "tests")

let fontPath = findFile("Arial.ttf")
doAssert(fontPath.len != 0, "Font file not found")

var rawData = readFile(fontPath)

var info: stbtt_fontinfo
if stbtt_InitFont(info, cast[ptr font_type](rawData.cstring), 0) == 0:
    doAssert(false, "Could not init font")

let b_w = 1024.cint # bitmap width
let b_h = 128.cint # bitmap height
let l_h = 64.cfloat # line height

# create a bitmap for the phrase
var bitmap = newSeq[byte](b_w * b_h)

# calculate font scaling
let scale = stbtt_ScaleForPixelHeight(info, l_h)

let word = toSeq("how are you? @ ¿ 1234567890 W".runes)

var x = 0

var ascent, descent, lineGap : cint
stbtt_GetFontVMetrics(info, ascent, descent, lineGap);

ascent = (ascent.cfloat * scale).cint
descent = (descent.cfloat * scale).cint

for i, c in word:
    # get bounding box for character (may be offset to account for chars that dip above or below the line
    var c_x1, c_y1, c_x2, c_y2 : cint
    stbtt_GetCodepointBitmapBox(info, c.cint, scale, scale, addr c_x1, addr c_y1, addr c_x2, addr c_y2)

    # compute y (different characters have different heights
    let y = ascent + c_y1

    # render character (stride and offset is important here)
    let byteOffset = x + (y  * b_w)
    stbtt_MakeCodepointBitmap(info, addr bitmap[byteOffset], c_x2 - c_x1, c_y2 - c_y1, b_w, scale, scale, word[i].cint)

    # how wide is this character
    var ax, lsb : cint
    stbtt_GetCodepointHMetrics(info, word[i].cint, ax, lsb)
    x += (ax.cfloat * scale).int

    # add kerning
    if i != word.len - 1:
        let kern = stbtt_GetCodepointKernAdvance(info, word[i].cint, word[i + 1].cint)
        x += (kern.cfloat * scale).int

proc saveBitmap(filename: string, data: openarray[byte], w, h: int) =
    var str = newString(data.len * 24)
    var i = 0
    for d in data:
        str[i] = chr(d)
        str[i + 1] = chr(d)
        str[i + 2] = chr(d)
        i += 3

    var s = newFileStream(filename, fmWrite)
    s.encodeBMP(str, w, h, 24)
    s.close()

let expectedPath = findFile("expected.bmp")
doAssert(expectedPath.len != 0, "expected.bmp not found")

saveBitmap("actual.bmp", bitmap, b_w, b_h)
doAssert(sameFileContent("actual.bmp", expectedPath))

# Test distance fields
import ../ttf/edtaa3func

let expectedDFPath = findFile("expected_df.bmp")
doAssert(expectedDFPath.len != 0, "expected.bmp not found")

make_distance_map(bitmap, b_w, b_h)
saveBitmap("actual_df.bmp", bitmap, b_w, b_h)
doAssert(sameFileContent("actual_df.bmp", expectedDFPath))

