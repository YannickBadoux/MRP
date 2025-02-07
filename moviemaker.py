## Made by Lucas Pouw, Yannick Badoux and Tim van der Vuurst for the 
## course "Simulations and Modeling in Astrophysics" '24-'25. 

import cv2
import os
import sys
import glob
import argparse
from tqdm import tqdm

# Move working directory one directory upwards (if running from the repo) to avoid saving anything in the github repo
if 'amuse-KL' in os.getcwd(): 
    os.chdir('..' ) 
print(f'Saving all files in: {os.getcwd()}\n')


def moviemaker(image_folder: str, video_name: str, fps: int =10, height: int =1080, width: int =1440, codec: str ='XVID') -> None:
    """ Uses cv2.VideoWriter to generate a video of plots generated by plotter.py. 

    Args:
        image_folder (str): Folder where images (that will act as frames) are stored.
        video_name (str): Name that is given to the made video.
        fps (int, optional): Frames played per second. Defaults to 10.
        height (int, optional): Height of the video in pixels. Defaults to 1080.
        width (int, optional): Width of the video in pixels. Defaults to 1440.
        codec (str, optional): codec to be passed to cv2.VideoWriter_fourcc. Defaults to 'XVID'. Recommended: avc1 or XVID
    """

    if not os.path.isdir(image_folder):
        sys.exit(f'{image_folder} is not an existing folder. Quitting.')

    print(f'Generating video {video_name} from {image_folder} with {fps} fps.')

    # Find all paths to images and sort them accordingly
    unsorted_images = glob.glob(image_folder+'/' + '*.png')
    image_numbers = [image_name.split('simulation')[-1].split('.')[0] for image_name in unsorted_images]
    image_numbers = list(map(int, image_numbers))
    _, sorted_images = zip(*sorted(zip(image_numbers, unsorted_images)))

    fourcc = cv2.VideoWriter_fourcc(*codec)
    video = cv2.VideoWriter(video_name, fourcc, fps, (width, height))

    for image_name in tqdm(sorted_images):
        img = cv2.imread(image_name)
        img = cv2.resize(img, (width, height))
        video.write(img)

    cv2.destroyAllWindows()
    video.release()


if __name__ == '__main__':
    # Create parser for modifying parameters in the terminal
    parser = argparse.ArgumentParser(description='Make a video from a folder of images.')
    parser.add_argument('--image_dir', type=str, required=True, help='Path to folder containing images. Required.')
    parser.add_argument('--video_name', type=str, required=True, help='Name of the video to be generated. Required.')
    parser.add_argument('--fps', default=10, type=int, help='Frames per second of the video.')
    parser.add_argument('--height', default=1080, type=int, help='Height of the video.')
    parser.add_argument('--width', default=1440, type=int, help='Width of the video.')
    parser.add_argument('--codec', default='XVID', type=str, help='Video codec used by OpenVC VideoWriter. (Recommended: avc1 or XVID)')
    args = parser.parse_args()

    #Pass along the args to the moviemaker function
    moviemaker(args.image_dir, args.video_name, args.fps, args.height, args.width, args.codec)
