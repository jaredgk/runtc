import unittest
import filecmp
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,'msh_est')))
from msh_from_vcf import getmsh
from aae_work import run_estimator

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


class estTest(unittest.TestCase):
    def test_sing_est(self):
        args = ['test_est_left_msh.txt','test_est_right_msh.txt','--n','50','--n0','10000','--mut','1e-8','--rec','1e-8','--round','3','--outfn','test_est.txt','--seed','1']
        run_estimator(args)
        self.assertTrue(filecmp.cmp('test_est.txt','test_sing_est.txt'))

    def test_sing_lowparams(self):
         args = ['test_est_left_msh.txt','test_est_right_msh.txt','--n','50','--n0','3000','--mut','1e-9','--rec','1e-9','--round','3','--outfn','test_est.txt','--seed','1']
         run_estimator(args)
         self.assertTrue(filecmp.cmp('test_est.txt','test_sing_lowparams_est.txt'))

    def test_alpha_est(self):
        args = ['test_alpha_est_left_msh.txt','test_alpha_est_right_msh.txt','--n','50','--n0','10000','--mut','1e-8','--rec','1e-8','--alledges','--round','3','--outfn','test_est.txt','--seed','1']
        run_estimator(args)
        self.assertTrue(filecmp.cmp('test_est.txt','test_est_alpha_est.txt'))

    def test_kall_est(self):
        args = ['test_kall_est_left_msh.txt','test_kall_est_right_msh.txt','--n','50','--n0','10000','--mut','1e-8','--rec','1e-8','--kmode','--round','3','--outfn','test_est.txt','--seed','1']
        run_estimator(args)
        self.assertTrue(filecmp.cmp('test_est.txt','test_est_kall_est.txt'))

#other tests: with genetic distances, estimator?

if __name__ == "__main__":
    unittest.main(verbosity=2)
        
