import unittest
import filecmp
import os
import sys
import subprocess

sys.path.insert(0, os.path.abspath(os.path.join(os.pardir,'msh_est')))
from msh_from_vcf import getmsh
from aae_work import run_estimator
import runtc

def callCCode(arglist):
    full_arglist = ['../msh_est/msh_vcf']+arglist
    print (full_arglist)
    s_call = subprocess.Popen(full_arglist,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    s_o,s_e = s_call.communicate()
    sys.stderr.write(s_e.decode())
    if s_call.returncode != 0:
        raise Exception("Error running c code\n%s"%(s_e.decode()))

class basicTest(unittest.TestCase):
    
    def test_alpha_withsgtn(self):
        args = ['--vcf','test_head.vcf','--out','test_head_ll.txt']
        getmsh(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_alpha_ws.txt'))
    
    def test_alpha_nosgtn(self):
        args = ['--vcf','test_head.vcf','--out','test_head_ll.txt','--exclude-singletons']
        getmsh(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_alpha_ns.txt'))

    def test_kall_withsgtn(self):
        args = ['--vcf','test_head.vcf','--k-all','--out','test_head_ll.txt']
        getmsh(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_kall_ws.txt'))

    def test_kall_nosgtn(self):
        args = ['--vcf','test_head.vcf','--k-all','--out','test_head_ll.txt','--exclude-singletons']
        getmsh(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_kall_ns.txt'))

    def test_pos(self):
        args = ['--vcf','test_head.vcf','--positions','test_pos.txt','--out','test_head_ll.txt']
        getmsh(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_pos.txt'))

    def test_sub(self):
        args = ['--vcf','test_head.vcf','--sub','test_sub.txt','--out','test_head_ll.txt']
        getmsh(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_sub.txt'))


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

    def test_twophase_est(self):
        args = ['test_est_left_msh.txt','test_est_right_msh.txt','--n','50','--n0','10000','--mut','1e-8','--rec','1e-8','--round','3','--outfn','test_est.txt','--seed','1','--twophase','.01','100']
        run_estimator(args)
        self.assertTrue(filecmp.cmp('test_est.txt','test_est_twophase_est.txt'))

class fullTest(unittest.TestCase):
    def test_basic(self):
        args = ['test_head.vcf','--map','test_map.txt','--mut','1e-8','--n0','10000','--outfn','test_est.txt','--outmsh','test_full_run','--k1']
        runtc.main(args)
        self.assertTrue(filecmp.cmp('test_est.txt','test_full_run.txt'))

    def test_cmsh(self):
        args = ['test_head.vcf','--map','test_map.txt','--mut','1e-8','--n0','10000','--outfn','test_est.txt','--outmsh','test_full_run','--c-msh','--k1']
        runtc.main(args)
        self.assertTrue(filecmp.cmp('test_est.txt','test_full_run_c.txt'))


class cppTest(unittest.TestCase):
    
    def test_alpha_withsgtn(self):
        args = ['--vcf','test_head.vcf','--out','test_head_ll.txt']
        callCCode(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_alpha_ws.txt'))
#other tests: with genetic distances, estimator?

    def test_alpha_nosgtn(self):
        args = ['--vcf','test_head.vcf','--exclude-singletons','--out','test_head_ll.txt']
        callCCode(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_alpha_ns.txt'))

    def test_kall_withsgtn(self):
        args = ['--vcf','test_head.vcf','--k-all','--out','test_head_ll.txt']
        callCCode(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_kall_ws.txt'))

    def test_kall_nosgtn(self):
        args = ['--vcf','test_head.vcf','--k-all','--exclude-singletons','--out','test_head_ll.txt']
        callCCode(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_kall_ns.txt'))

    def test_pos(self):
        args = ['--vcf','test_head.vcf','--positions','test_pos.txt','--out','test_head_ll.txt']
        callCCode(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_pos.txt'))

    def test_sub(self):
        args = ['--vcf','test_head.vcf','--sub','test_sub.txt','--out','test_head_ll.txt']
        callCCode(args)
        self.assertTrue(filecmp.cmp('test_head_ll.txt','test_head_sub.txt'))

#other tests: with genetic distances, estimator?

if __name__ == "__main__":
    unittest.main(verbosity=2)
        
