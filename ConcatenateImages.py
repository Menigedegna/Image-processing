# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 05:47:58 2017

@author: Pheonix
"""

from PIL import Image
import os

path=r"C:\Users\Pheonix\Desktop\Test.tif.tif"
result = Image.new("RGB", (800, 800))

files = [path, path,path,path]

for index, file in enumerate(files):
    path = os.path.expanduser(file)
    img = Image.open(path)
    img.thumbnail((400, 400), Image.ANTIALIAS)
    x = index // 2 * 400
    y = index % 2 * 400
    w, h = img.size
    print('pos {0},{1} size {2},{3}'.format(x, y, w, h))
    result.paste(img, (x, y, x + w, y + h))
result.save(os.path.expanduser(r"C:\Users\Pheonix\Desktop\Result.tif"))