import pytest

if __name__ == '__main__':
    # It's expected that we are being run from the test subdirectory,
    # and so all the test modules and the mesh directory is under
    # 'validation/'.
    pytest.main(args=[__file__,'-s','--all-modules','-v','validation_rd'])
    pytest.main(args=[__file__,'-s','--all-modules','-v','validation_cp'])
    pytest.main(args=[__file__,'-s','--all-modules','-v','validation_efield'])

