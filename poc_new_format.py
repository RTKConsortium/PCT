# Dummy example of usage of the new PCT format,
# which defines a dataset that only contains upstream position of the protons.
import itk
import numpy as np

TEST_FILE = "/tmp/test.mhd"

# Create a list-mode dataset that uses the new format
data = itk.image_from_array(np.random.rand(10, 3))
data["UpstreamPositionU"] = 0
data["UpstreamPositionV"] = 1
data["UpstreamPositionW"] = 2

itk.imwrite(data, TEST_FILE)


# Read from a file that uses the new format
class ProtonPairs:

    pct_fields = [
        "UpstreamPositionU",
        "UpstreamPositionV",
        "UpstreamPositionW",
        "DownstreamPositionU",
        "DownstreamPositionV",
        "DownstreamPositionW",
    ]

    def __init__(self, file_path):
        im = itk.imread(file_path)

        self.data = itk.array_from_image(im)
        self.metadata = dict(im)

    def __getitem__(self, idx):
        # __getitem__ is enough to implement the iterator protocol
        data_idx = self.data[idx]
        item = {}
        for pct_field in ProtonPairs.pct_fields:
            try:
                item[pct_field] = data_idx[int(self.metadata[pct_field])]
            except KeyError:
                pass
        return item


proton_pairs = ProtonPairs(TEST_FILE)

for i, p in enumerate(proton_pairs):
    print("Proton pair", i, ": ", p)

p_iter = iter(proton_pairs)
print(next(p_iter))
print(next(p_iter))
print(next(p_iter))
