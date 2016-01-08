import sequtils, math

# Original is taken from freetype-gl

{.push checks: off, stackTrace: off.}

# Compute the local gradient at edge pixels using convolution filters.
# The gradient is computed only at edge pixels. At other places in the
# image, it is never used, and it's mostly zero anyway.
proc computegradient(img: openarray[cdouble], w, h: cint, gx, gy: var openarray[cdouble]) =
    const SQRT2 = cdouble(1.4142136)
    for i in 1 ..< h - 1: # Avoid edges where the kernels would spill over
        for j in 1 ..< w - 1:
            let k = i * w + j
            if (img[k] > 0) and (img[k] < 1): # Compute gradient for edge pixels only
                gx[k] = -img[k-w-1] - SQRT2*img[k-1] - img[k+w-1] + img[k-w+1] + SQRT2*img[k+1] + img[k+w+1]
                gy[k] = -img[k-w-1] - SQRT2*img[k-w] - img[k+w-1] + img[k-w+1] + SQRT2*img[k+w] + img[k+w+1]
                var glength = gx[k]*gx[k] + gy[k]*gy[k]
                if glength > 0: # Avoid division by zero
                    glength = sqrt(glength)
                    gx[k] /= glength
                    gy[k] /= glength
    # TODO: Compute reasonable values for gx, gy also around the image edges.
    # (These are zero now, which reduces the accuracy for a 1-pixel wide region
    # around the image edge.) 2x2 kernels would be suitable for this.

# A somewhat tricky function to approximate the distance to an edge in a
# certain pixel, with consideration to either the local gradient (gx,gy)
# or the direction to the pixel (dx,dy) and the pixel greyscale value a.
# The latter alternative, using (dx,dy), is the metric used by edtaa2().
# Using a local estimate of the edge gradient (gx,gy) yields much better
# accuracy at and near edges, and reduces the error even at distant pixels
# provided that the gradient direction is accurately estimated.
proc edgedf(gx, gy, a: cdouble): cdouble =
    if (gx == 0) or (gy == 0): # Either A) gx or gv are zero, or B) both
        result = 0.5 - a  # Linear approximation is A) correct or B) a fair guess
    else:
        let glength = sqrt(gx*gx + gy*gy)
        var ggx = gx
        var ggy = gy
        if glength > 0:
            ggx /= glength
            ggy /= glength
        # Everything is symmetric wrt sign and transposition,
        # so move to first octant (ggx>=0, ggy>=0, ggx>=ggy) to
        # avoid handling all possible edge directions.

        ggx = abs(ggx)
        ggy = abs(ggy)
        if ggx < ggy:
            swap(ggx, ggy)
        let a1 = 0.5 * ggy / ggx
        if a < a1: # 0 <= a < a1
            result = 0.5 * (ggx + ggy) - sqrt(2.0 * ggx * ggy * a)
        elif a < 1 - a1: # a1 <= a <= 1-a1
            result = (0.5 - a) * ggx
        else: # 1-a1 < a <= 1
            result = -0.5 * (ggx + ggy) + sqrt(2.0 * ggx * ggy * (1.0 - a))

proc distaa3(img, gximg, gyimg: openarray[cdouble], w, c, xc, yc, xi, yi: cint): cdouble =
    let closest = c-xc-yc*w # Index to the edge pixel pointed to from c
    var a = img[closest]    # Grayscale value at the edge pixel
    let gx = gximg[closest] # X gradient component at the edge pixel
    let gy = gyimg[closest] # Y gradient component at the edge pixel

    if(a > 1.0): a = 1.0
    if(a < 0.0): a = 0.0 # Clip grayscale values outside the range [0,1]
    if(a == 0.0): return 1000000.0 # Not an object pixel, return "very far" ("don't know yet")

    let dx = xi.cdouble
    let dy = yi.cdouble
    let di = sqrt(dx*dx + dy*dy) # Length of integer vector, like a traditional EDT
    if(di==0): # Use local gradient only at edges
        # Estimate based on local gradient only
        result = edgedf(gx, gy, a)
    else:
        # Estimate gradient based on direction to edge (accurate for large di)
        result = edgedf(dx, dy, a)
    result += di # Same metric as edtaa2, except at edges (where di=0)

proc edtaa3(img, gx, gy: openarray[cdouble], w, h: cint, distx, disty: var openarray[int16], dist: var openarray[cdouble]) =
    var c : int
    var olddist, newdist: cdouble
    var cdistx, cdisty: int16
    var newdistx, newdisty: int16

    const epsilon = 1e-3

    # Shorthand template: add ubiquitous parameters dist, gx, gy, img and w and call distaa3()
    template DISTAA(c, xc, yc, xi, yi): cdouble = distaa3(img, gx, gy, w, c.cint, xc.cint, yc.cint, xi.cint, yi.cint)

    # Initialize index offsets for the current image width
    let offset_u = -w
    let offset_ur = -w + 1
    let offset_r = 1
    let offset_rd = w + 1
    let offset_d = w
    let offset_dl = w-1
    let offset_l = -1
    let offset_lu = -w - 1

    # Initialize the distance images
    let sz = w * h
    var i = 0
    while i < sz:
        distx[i] = 0 # At first, all pixels point to
        disty[i] = 0 # themselves as the closest known.
        if img[i] <= 0:
            dist[i]= 1000000 # Big value, means "not set yet"
        elif img[i] < 1:
            dist[i] = edgedf(gx[i], gy[i], img[i]) # Gradient-assisted estimate
        else:
            dist[i]= 0 # Inside the object
        inc i

    var changed = true

    # Perform the transformation
    while changed: # Sweep until no more updates are made
        changed = false
        # Scan rows, except first row
        var y = 1
        while y < h:
            # move index to leftmost pixel of current row
            var i = y * w
            # scan right, propagate distances from above & left

            # Leftmost pixel is special, has no left neighbors
            olddist = dist[i]
            if olddist > 0: # If non-zero distance or not set yet
                var c = i + offset_u # Index of candidate for testing
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx
                newdisty = cdisty+1
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    olddist=newdist
                    changed = true

                c = i+offset_ur
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx-1
                newdisty = cdisty+1
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    changed = true
            inc i

            var x = 1
            # Middle pixels have all neighbors
            while x < w - 1:
                defer:
                    inc x
                    inc i
                olddist = dist[i]
                if olddist <= 0: continue # No need to update further
                c = i+offset_l
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx+1
                newdisty = cdisty
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    olddist=newdist
                    changed = true

                c = i+offset_lu
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx+1
                newdisty = cdisty+1
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if(newdist < olddist-epsilon):
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    olddist=newdist
                    changed = true

                c = i+offset_u
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx
                newdisty = cdisty+1
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    olddist=newdist
                    changed = true

                c = i+offset_ur
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx-1
                newdisty = cdisty+1
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    changed = true

            # Rightmost pixel of row is special, has no right neighbors
            olddist = dist[i]
            if olddist > 0: # If not already zero distance
                c = i+offset_l
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx+1
                newdisty = cdisty
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    olddist=newdist
                    changed = true

                c = i+offset_lu
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx+1
                newdisty = cdisty+1
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    olddist=newdist
                    changed = true

                c = i+offset_u
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx
                newdisty = cdisty+1
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    changed = true

            # Move index to second rightmost pixel of current row.
            # Rightmost pixel is skipped, it has no right neighbor.
            i = y*w + w-2

            # scan left, propagate distance from right
            x = w-2
            while x >= 0:
                defer:
                    dec x
                    dec i
                olddist = dist[i]
                if olddist <= 0: continue # Already zero distance

                let c = i+offset_r
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx-1
                newdisty = cdisty
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    changed = true
            inc y

        # Scan rows in reverse order, except last row
        y = h - 2
        while y >= 0:
            defer: dec y
            # move index to rightmost pixel of current row
            var i = y*w + w-1

            # Scan left, propagate distances from below & right

            # Rightmost pixel is special, has no right neighbors
            olddist = dist[i]
            if olddist > 0: # If not already zero distance
                c = i+offset_d
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx
                newdisty = cdisty-1
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    olddist=newdist
                    changed = true

                c = i+offset_dl
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx+1
                newdisty = cdisty-1
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    changed = true
            dec i

            # Middle pixels have all neighbors
            var x = w - 2
            while x > 0:
                defer:
                    dec x
                    dec i
                olddist = dist[i]
                if olddist <= 0: continue # Already zero distance

                c = i+offset_r
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx-1
                newdisty = cdisty
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    olddist=newdist
                    changed = true

                c = i+offset_rd
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx-1
                newdisty = cdisty-1
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    olddist=newdist
                    changed = true

                c = i+offset_d
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx
                newdisty = cdisty-1
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    olddist=newdist
                    changed = true

                c = i+offset_dl
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx+1
                newdisty = cdisty-1
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    changed = true

            # Leftmost pixel is special, has no left neighbors
            olddist = dist[i]
            if olddist > 0: # If not already zero distance
                c = i+offset_r
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx-1
                newdisty = cdisty
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    olddist=newdist
                    changed = true

                c = i+offset_rd
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx-1
                newdisty = cdisty-1
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    olddist=newdist
                    changed = true

                c = i+offset_d
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx
                newdisty = cdisty-1
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    changed = true


            # Move index to second leftmost pixel of current row.
            # Leftmost pixel is skipped, it has no left neighbor.
            i = y*w + 1
            x = 1
            while x < w:
                defer:
                    inc x
                    inc i
                # scan right, propagate distance from left */
                olddist = dist[i]
                if olddist <= 0: continue # Already zero distance

                c = i+offset_l
                cdistx = distx[c]
                cdisty = disty[c]
                newdistx = cdistx+1
                newdisty = cdisty
                newdist = DISTAA(c, cdistx, cdisty, newdistx, newdisty)
                if newdist < olddist-epsilon:
                    distx[i]=newdistx
                    disty[i]=newdisty
                    dist[i]=newdist
                    changed = true

when defined(js):
    proc newSeqDouble(sz: int): seq[cdouble] {.importc: "new Float32Array".}
    proc newSeqInt16(sz: int): seq[int16] {.importc: "new Int16Array".}
else:
    template newSeqDouble(sz: int): seq[cdouble] = newSeq[cdouble](sz)
    template newSeqInt16(sz: int): seq[int16] = newSeq[int16](sz)

proc make_distance_map*(data: var openarray[cdouble], width, height : cuint) =
    let sz = (width * height).int
    var xdist = newSeqInt16(sz)
    var ydist = newSeqInt16(sz)
    var gx = newSeqDouble(sz)
    var gy = newSeqDouble(sz)
    var outside = newSeqDouble(sz)
    var inside = newSeqDouble(sz)

    var img_min = 255.cdouble
    var img_max = -255.cdouble

    # Convert img into double (data)
    var i = 0
    while i < sz:
        let v = data[i]
        if v > img_max: img_max = v
        if v < img_min: img_min = v
        inc i

    # Rescale image levels between 0 and 1
    i = 0
    while i < sz:
        data[i] = (data[i] - img_min) / img_max
        inc i

    # Compute outside = edtaa3(bitmap); % Transform background (0's)
    computegradient(data, width.cint, height.cint, gx, gy)
    edtaa3(data, gx, gy, width.cint, height.cint, xdist, ydist, outside)

    i = 0
    while i < sz:
        if outside[i] < 0:
            outside[i] = 0
        inc i

    # Compute inside = edtaa3(1-bitmap); % Transform foreground (1's)
    gx.applyIt(0)
    gy.applyIt(0)

    i = 0
    while i < sz:
        data[i] = 1 - data[i]
        inc i

    computegradient(data, width.cint, height.cint, gx, gy)
    edtaa3(data, gx, gy, width.cint, height.cint, xdist, ydist, inside)

    i = 0
    while i < sz:
        if inside[i] < 0:
            inside[i] = 0
        inc i

    # distmap = outside - inside; % Bipolar distance field
    i = 0
    while i < sz:
        outside[i] -= inside[i]
        outside[i] = 128+outside[i]*16
        if outside[i] < 0: outside[i] = 0
        if outside[i] > 255: outside[i] = 255
        data[i] = outside[i]
        inc i

proc make_distance_map*(img: var openarray[byte], width, height : cuint) =
    let sz = (width * height).int
    var data = newSeqDouble(sz)

    # Convert img into double (data)
    var i = 0
    while i < sz:
        data[i] = img[i].cdouble
        inc i

    make_distance_map(data, width, height)

    # Convert back
    i = 0
    while i < sz:
        img[i] = (255.byte - data[i].byte)
        inc i

{.pop.}
