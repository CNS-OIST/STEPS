import nose

if __name__ == '__main__':
    # It's expected that we are being run from the test subdirectory,
    # and so all the test modules and the mesh directory is under
    # 'validation/'.

    nose.run(argv=[__file__,'-s','--all-modules','-v','validation_cp'])
    nose.run(argv=[__file__,'-s','--all-modules','-v','validation_efield'])
    nose.run(argv=[__file__,'-s','--all-modules','-v','validation_rd'])



