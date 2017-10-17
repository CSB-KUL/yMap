from __future__ import print_function, absolute_import, division, unicode_literals
import unittest
import tempfile
import shutil
import os
import time
import filecmp
from ymap.ymap import data, ymap_proteins, ymap_genes, web, YGtPM

ref_dir = os.path.abspath('test_files')


class Timer:
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start


class YGtPMTest(unittest.TestCase):
    """
    YMap tests
    """

    def setUp(self):
        # Create a temporary directory
        # self.ref_dir = os.path.abspath('test_files')
        self.test_dir = tempfile.mkdtemp()
        os.chdir(self.test_dir)
        self.c = YGtPM()
        self.store_ref = False
        # at various places out put is compared to stored reference files. put this to true
        # to regenerate all reference files

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)

    def test_class_creation(self):
        """
        testing the creation of YGtPM
        """
        self.assertIsInstance(self.c, YGtPM)

    def test_data(self):
        """
        testing the steps in data()
        """
        self.assertIsInstance(self.c, YGtPM)
        try:  # the actual download can only be performed a limit number of times per hour
            raise RuntimeError  # this is to prevent connecting the database during development comment for real tesing
            self.c.pTMdata()
            if self.store_ref:
                shutil.copy('uniprot_mod_raw.txt', ref_dir)  # saving the reference data
        except:
            shutil.copy(os.path.join(ref_dir, 'uniprot_mod_raw.txt'), '.')
        self.assertTrue(os.path.isfile('uniprot_mod_raw.txt'))

        self.c.clean('uniprot_mod_raw.txt')  # produces PTMs.txt
        self.assertTrue(os.path.isfile('PTMs.txt'))
        if self.store_ref:
            shutil.copy('PTMs.txt', ref_dir)  # saving the reference data
        self.assertTrue(filecmp.cmp('PTMs.txt', os.path.join(ref_dir, 'PTMs.txt')))

        self.c.iD()
        self.assertTrue(os.path.isfile('yeastID.txt'))
        if self.store_ref:
            shutil.copy('yeastID.txt', ref_dir)  # saving the reference data
        self.assertTrue(filecmp.cmp('yeastID.txt', os.path.join(ref_dir, 'yeastID.txt')))

        self.c.pmap('yeastID.txt', 'PTMs.txt')
        self.assertTrue(os.path.isfile('PTM_id_file.txt'))
        if self.store_ref:
            shutil.copy('PTM_id_file.txt', ref_dir)  # saving the reference data
        self.assertTrue(filecmp.cmp('PTM_id_file.txt', os.path.join(ref_dir, 'PTM_id_file.txt')))
        self.assertTrue(False)

    def test_ymap_proteins(self):
        self.assertTrue(True)

    def test_ymap_genes(self):
        self.assertTrue(True)

    def test_web(self):
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()
