import sys
import os
from PIL import Image

def change_size(path_in, path_out):
    img = Image.open(path_in)
    pal = img.getpalette()
    width, height = img.size

    if width >= height:
        newsize = [width, width]
        buf_height = int((width-height)/2)
        buf_width = 0
    else:
        newsize = [height, height]
        buf_height = 0
        buf_width = int((height-width)/2)
    result = Image.new('LA', newsize)

    # 压缩
    im = img.load()
    print(im[0,0])
    res = result.load()
    for x in range(width):
        for y in range(height):
            t = im[x, y][-1]
            color = im[x, y][0]
            res[x+buf_width, y+buf_height] = (color, t)

    result.resize((100, 100), Image.ANTIALIAS).save(path_out)

for file in os.listdir('pics_before'):
    print(file)
    path_in = 'pics_before/'+file
    path_out = 'pics_after/'+file
    change_size(path_in, path_out)