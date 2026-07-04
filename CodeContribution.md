# Code contribution

## Coding style

PCT is based on RTK and ITK and aims at following their coding conventions. Any developer should follow these conventions when submitting new code or contributions to the existing one. We strongly recommend reading thoroughly [ITK's style guide](http://www.itk.org/Wiki/ITK/Coding_Style_Guide) and [RTK coding conventions](https://docs.openrtk.org/en/latest/CodeContribution.html).

## Testing

This section describes how to add/edit datasets for testing purposes for PCT.

### Test data

Datasets are not stored in the GIT repository for efficiency and also to avoid having large history due to binary files. Instead, the files are stored on a [Girder](http://data.kitware.com) instance. Here's the recipe to add new datasets:

1.  Register/Login to Girder hosted at Kitware: [http://data.kitware.com](http://data.kitware.com)
2.  Locate the RTK collection: [https://data.kitware.com/#collection/5a7706878d777f0649e04776](https://data.kitware.com/#collection/5a7706878d777f0649e04776)
3.  Upload the new datasets in the folder named`PCTTestingData`. If you do not have the necessary privileges, please email the [RTK mailing list](https://www.creatis.insa-lyon.fr/mailman/listinfo/rtk-users)
4.  In the GIT repository, add a file in `test/Baseline` or `test/Input` with the exact filename as the original file **but with the .md5 extension**. Inside that file put the md5sum of the file as found on Girder.

### C++ test

Add a call to the `itk_add_test` macro in `test/CMakeLists.txt` and specify the datasets you want CTest to download by appending the data to `DATA{}`. For example:

```cmake
itk_add_test(NAME pctapppaircutstest
    COMMAND itkTestDriver
      --compare
      DATA{Baseline/baseline_paircuts.mhd,baseline_paircuts.raw}
      ${ITK_TEST_OUTPUT_DIR}/paircuts.mhd
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/pctpaircuts
      -i DATA{Input/input_paircuts.mhd,input_paircuts.raw} -o ${ITK_TEST_OUTPUT_DIR}/paircuts.mhd --roiradius 200 --minwepl 70 --maxwepl 180 --fluence .5 --seed 1234
  )
```

### Python test

PCT uses [`pytest`](https://docs.pytest.org/en/stable/) for the Python tests and follows its conventions. To add a Python test, create a Python file in the `test` folder. The name of the file should end in `_test.py` in order to be automatically picked up by `pytest`. Then create a test function whose name starts with `test_`. This function will be automatically executed when running `pytest`. For example:

```python
def test_python_wrapping_instantiation():
    imageType = itk.Image[itk.F, 3]
    pct.FDKDDWeightProjectionFilter[imageType].New()
    pct.FDKDDConeBeamReconstructionFilter[imageType].New()
    pct.FDKDDBackProjectionImageFilter[imageType, imageType].New()
    pct.ProtonPairsToDistanceDrivenProjection[imageType, imageType].New()
    pct.HoleFillingImageFilter[imageType].New()
```
