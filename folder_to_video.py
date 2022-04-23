import cv2
import glob

img_array = []
for filename in sorted(glob.glob("animate/*.png"), key=lambda x: int((x.split("\\")[-1]).split(".")[0])):
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width, height)
    img_array.append(img)

out = cv2.VideoWriter('project.avi', cv2.VideoWriter_fourcc(*'DIVX'), 15, size)

for i in range(len(img_array)):
    out.write(img_array[i])
out.release()