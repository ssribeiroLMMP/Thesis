import cv2
import numpy as np
import os

from os.path import isfile, join

def convert_frames_to_video(pathIn,pathOut,fps,textNum):
    frame_array = []
    files = [f for f in os.listdir(pathIn) if isfile(join(pathIn, f))]

    #for sorting the file names properly
    files.sort(key = lambda x: int(x[textNum:-4]))

    for i in range(len(files)):
        filename=pathIn + files[i]
        #reading each files
        img = cv2.imread(filename)
        height, width, layers = img.shape
        size = (width,height)
        print(filename)
        #inserting the frames into an image array
        frame_array.append(img)

    out = cv2.VideoWriter(pathOut,cv2.VideoWriter_fourcc(*'DIVX'), fps, size)

    for i in range(len(frame_array)):
        # writing to a image array
        out.write(frame_array[i])
    out.release()

def main():
    pathIn= './/Images//'
    #'.//Pendant - Air_DI Water, 72.0mN_m - SM500microL, N22, 5microL, 20ºC, NB0.14, 3600s - 26.8.2019//' # Heloisa - ADIW.avi
    #'.//Rising - DI water_hexadecane, 44.2mN_m - SM500microL, N18, 32microL, 25ºC, NB0.16, 61200s - 27.8.2019//'# Heloisa - DIWH.avi
    #'.//Rising - Sea water_hexadecane, 26.9mN_m - SM500microL, N18, 22microL, 25ºC, NB0.14, 25h - 28.8.2019//'  # Heloisa - SWH.avi 
    pathOut = 'FEniCs_HeleShaw2D_VR=0_8.avi'
    fps = 12.0
    textNum = len('Concentration.')
    convert_frames_to_video(pathIn, pathOut, fps, textNum)

if __name__=="__main__":
    main()
