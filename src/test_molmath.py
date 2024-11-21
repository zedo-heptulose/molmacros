import unittest
import molmath as mm
import numpy as np

abtests =  [
        [np.array([1,0,0]),np.array([0,1,0]),np.array([0,0,1]),np.pi/2],
        [np.array([0,1,0]),np.array([0,1,0]),np.array([0,0,1]),0],
        [np.array([0,-1,0]),np.array([0,1,0]),np.array([0,0,1]),np.pi],
        [np.array([0.5,np.sqrt(3)/2,0]),np.array([0,1,0]),np.array([0,0,1]),np.pi/6],
        [np.array([np.sqrt(3)/2,-0.5,0]),np.array([0,1,0]),np.array([0,0,1]),2*np.pi/3],
        
        [np.array([0,1,0]),np.array([1,0,0]),np.array([0,0,1]),2*np.pi-np.pi/2],
        [np.array([0,1,0]),np.array([0,1,0]),np.array([0,0,1]),2*np.pi-0],
        [np.array([0,1,0]),np.array([0,-1,0]),np.array([0,0,1]),2*np.pi-np.pi],
        [np.array([0,1,0]),np.array([0.5,np.sqrt(3)/2,0]),np.array([0,0,1]),2*np.pi -np.pi/6],
        [np.array([0,1,0]),np.array([np.sqrt(3)/2,-0.5,0]),np.array([0,0,1]),2*np.pi-2*np.pi/3]
    ]

class TestMolMath(unittest.TestCase):
    def test_angle_between(self):
        print('\n\n//////////test_angle_between///////////\n\n')
        for test in abtests:
            print(f'test: {test}')
            self.assertAlmostEqual(mm.angle_about_axis(test[0],test[1],test[2]),test[3]%(2*np.pi))
                  
    def test_rotation_about_axis(self):
        print('\n\n//////////test_rotation_about_axis///////////\n\n')
        for test in abtests:
            print(f'test: {test}')
            v1 = test[0]
            v2 = test[1]
            axis = test[2]
            angle = test[3]
            rotation_matrix = mm.rotation_about_axis(axis, angle)
            print(f'rotation_matrix: \n{rotation_matrix}\n\n')
            self.assertTrue(np.allclose(rotation_matrix @ v1, v2,atol=1e-2))
                  
    def test_align_matrix(self):
        print('\n\n//////////test_align_matrix///////////\n\n')
        for test in abtests:
            print(f'test: {test}')
            v1 = test[0]
            v2 = test[1]
            axis = test[2]
            angle = test[3]
            align_matrix = mm.align_matrix(v1, v2)
            print(f'align_matrix: \n{align_matrix}\n\n')
            self.assertTrue(np.allclose(align_matrix @ v2, v1 ,atol=1e-2))
            # rotation_matrix = mm  .rotation_about_axis(axis, angle)
            # self.assertTrue(np.allclose(align_matrix @ v1, rotation_matrix @ v1,atol=1e-2))
            # self.assertTrue(np.allclose(v1, v3, atol=1e-2))

if __name__ == '__main__':
    unittest.main()
