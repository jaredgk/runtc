import unittest
import filecmp
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,'msh_est')))
from msh_from_vcf import getmsh

class basicTest(unittest.TestCase):
    
    def test_alpha_withsgtn(self):
        args = ['--vcf','test_head.vcf','--out','test_head_ll.txt']
        getmsh(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_alpha_ws.txt'))
    
    def test_alpha_nosgtn(self):
        args = ['--vcf','test_head.vcf','--out','test_head_ll.txt','--exclude-singletons']
        getmsh(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_alpha_ns.txt'))

    def test_singleton_withsgtn(self):
        args = ['--vcf','test_head.vcf','--singleton','--out','test_head_ll.txt']
        getmsh(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_sing_ws.txt'))

    def test_singleton_nosgtn(self):
        args = ['--vcf','test_head.vcf','--singleton','--out','test_head_ll.txt','--exclude-singletons']
        getmsh(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_sing_ns.txt'))

    def test_kall_withsgtn(self):
        args = ['--vcf','test_head.vcf','--k-all','--out','test_head_ll.txt']
        getmsh(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_kall_ws.txt'))

    def test_kall_nosgtn(self):
        args = ['--vcf','test_head.vcf','--k-all','--out','test_head_ll.txt','--exclude-singletons']
        getmsh(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_kall_ns.txt'))

#other tests: with genetic distances, estimator?

if __name__ == "__main__":
    unittest.main(verbosity=2)
        
