import os
from typing import Mapping, Tuple, Union
import numpy as np
from starfish.types import Axes, Coordinates, Features, Number
from starfish import Codebook
from starfish.experiment.builder import FetchedTile, TileFetcher
from slicedimage import ImageFormat
from skimage.io import imread
from starfish.experiment.builder import write_experiment_json

# imaging
channels = ["DAPI", "AF488", "Cy3", "Cy5", "AF750", "Atto425"]
pixelscale = 0.1625
zstep = 0.41

# codebook in csv (column 1: code, column2-x: color code in each round)
# AVOID comma in gene names (column 1)
codebook_csv = r"D:\HCA_07_DG\human_codebook.csv"
DO_decorators = ["Atto425" "AF488" "Cy3" "Cy5"]

# tile position
# TODO: read czi metadata using bio-format
tilepos_xy_csv = r"D:\HCA_07_DG\human_tile_coordinates.csv"

# not sure why this is so important and can't be read from image itself...
tilesz = 2048
SHAPE = tilesz, tilesz

# file tiff location
input_dir = r"D:\HCA_07_DG\human_test_1FOV\3D"
output_dir = r"D:\HCA_07_DG\starfish_190208\3D"


# hakuna matata
if not os.path.isdir(output_dir):
    try:
        os.mkdir(output_dir)
    except:
        os.makedirs(output_dir)

# TODO: register, re-slice and format images


def make_codebook_json(codebook_csv):
    """ convert color code matrix in csv to json format"""
    codebook_array = []
    with open(codebook_csv, "r") as f:
        for line in f:
            line = line.rstrip('\n').split(',')
            codewords = []
            for r, colorcode in enumerate(line[1:]):
                codewords.append({Axes.ROUND.value: r,
                                  Axes.CH.value: int(colorcode)-1,
                                  Features.CODE_VALUE:1})
            codebook_array.append({Features.CODEWORD:codewords, Features.TARGET: line[0]})
    codebook = Codebook.from_code_array(codebook_array)
    codebook_json_filename = "codebook.json"
    codebook.to_json(os.path.join(output_dir, codebook_json_filename))


def add_codebook(experiment_json_doc):
    experiment_json_doc['codebook'] = "codebook.json"

    return experiment_json_doc


def get_tilepos(tilepos_xy_csv):
    tilexy = np.ndarray((0,2))
    with open(tilepos_xy_csv, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split(',')
            tilexy = np.vstack([tilexy, [np.double(line[0]), np.double(line[1])]])
    return tilexy


class ISSTile3D(FetchedTile):
    def __init__(self, file_path, fov, z):
        self.file_path = file_path
        self.fov = fov
        self.z = z

    @property
    def shape(self) -> Tuple[int, ...]:
        return SHAPE

    @property
    def coordinates(self) -> Mapping[Union[str, Coordinates], Union[Number, Tuple[Number, Number]]]:
        # FIXME: (dganguli) please provide proper coordinates here.
        return {
            Coordinates.X: (tilexy[self.fov, 0]*pixelscale, (tilexy[self.fov, 0] + tilesz)*pixelscale),
            Coordinates.Y: (tilexy[self.fov, 1]*pixelscale, (tilexy[self.fov, 1] + tilesz)*pixelscale),
            Coordinates.Z: (self.z*zstep, (self.z+1)*zstep),
        }

    def tile_data(self) -> np.ndarray:
        return imread(self.file_path)


class ISS3DPrimaryTileFetcher(TileFetcher):
    def __init__(self, input_dir):
        self.input_dir = input_dir

    def get_tile(self, fov: int, r: int, ch: int, z: int) -> FetchedTile:
        # return ISSTile3D(os.path.join(self.input_dir, "r{}_c{}_z{}_t{}.tif".format(r+1, ch+1, z+1, fov+1)), fov, z)
        return ISSTile3D(os.path.join(self.input_dir, "r{}_c{}_z{}_t{}.tif".format(r+1, ch+1, z+1, 340+1)), 340, z)


class ISS3DAuxTileFetcher(TileFetcher):
    def __init__(self, path):
        self.path = path

    def get_tile(self, fov: int, r: int, ch: int, z: int) -> FetchedTile:
        return ISSTile3D(os.path.join(self.path, "r{}_c0_z{}_t{}.tif".format(r+1, z+1, 340+1)), 340, z)


tilexy = get_tilepos(tilepos_xy_csv)
make_codebook_json(codebook_csv)


write_experiment_json(
    path=output_dir, fov_count=1, tile_format=ImageFormat.TIFF,
    primary_image_dimensions={
        Axes.ROUND: 4,
        Axes.CH: 4,
        Axes.ZPLANE: 18,
    },
    aux_name_to_dimensions={
        'nuclei': {
            Axes.ROUND: 1,
            Axes.CH: 1,
            Axes.ZPLANE: 18,
        },
        # 'dots': {
        #     Axes.ROUND: 1,
        #     Axes.CH: 1,
        #     Axes.ZPLANE: 1,
        # },
    },
    primary_tile_fetcher=ISS3DPrimaryTileFetcher(input_dir),
    aux_tile_fetcher={
        'nuclei': ISS3DAuxTileFetcher(input_dir),
    },
    postprocess_func=add_codebook,
    default_shape=SHAPE
)
