import numpy as np
import json
from LaPy.lapy.TriaIO import import_off
import os
from sklearn.cluster import KMeans
from itertools import permutations, combinations
import pandas as pd

class ProcrustesComp:
    def __init__(self):
        self.X = np.identity(3)
        self.Y = np.identity(3)
        self.Z = np.identity(3)

    def MSE(self, A, B):
        if A.shape == B.T.shape:
            N = A.shape[0] * A.shape[1]
            diff = np.subtract(A, B.T)
            mse = np.sum(np.power(diff, 2))/N
            return mse
        elif A.shape[0] < B.shape[0]:
            B_kmm = KMeans(n_clusters=A.shape[0], random_state=0).fit(B)
            N = A.shape[0] * A.shape[1]
            diff = np.subtract(A, B_kmm.cluster_centers_.T)
            mse = np.sum(np.power(diff, 2))/N
            return mse

    def takeToPreshape(self, lm_dict):
        #keys = list(lm_dict.keys())
        #landmark_array_full = np.array([lm_dict[keys[i]] for i in range(len(keys))])
        landmark_array_nn = lm_dict[~np.isnan(lm_dict).any(axis=1), :]
        preshape = (landmark_array_nn - np.mean(landmark_array_nn,axis=0))/np.linalg.norm(landmark_array_nn, ord='fro')
        return preshape

    def setRotationMatrices(self, axis, theta):
        if axis =='x':
            self.X = np.array([[1,0,0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
        elif axis=='y':
            self.Y = np.array([[np.cos(theta), 0, -np.sin(theta)], [0,1,0], [np.sin(theta), 0, np.cos(theta)]])
        else:
            self.Z = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0,0,1]])
    def rotateShape(self, B):
        rot = self.X @ self.Y @ self.Z 
        return rot @ B.T
        
    def compare(self, lm_dict_A, lm_dict_B):
        psA = self.takeToPreshape(lm_dict_A)
        psB = self.takeToPreshape(lm_dict_B)
        if psA.shape[0] != psB.shape[0]:
            # PRUNE
            mse = self.prune(psA, psB)
        else:
            # NO PRUNE, JUST SEARCH
            mse = self.search(psA, psB)
        return mse

    def findBestMatch(self, poses_dict, B):
        res = []
        for p in poses_dict.keys():
            print(f'Evaluating against pose {p}')
            A = poses_dict[p]
            mse = self.compare(A, B)
            res.append(mse)
        return res

    def search(self, dict_landmarks, target_landmarks):
        xy_theta = [(np.pi/32)*i for i in range(-6, 7)]
        z_theta = [i*np.pi/32 for i in range(65)]
        mseMin = 1e10
        best = []
        for x_thet in xy_theta:
            self.setRotationMatrices('x', x_thet)
            for y_thet in xy_theta:
                self.setRotationMatrices('y', y_thet)
                for z_thet in z_theta:
                    self.setRotationMatrices('z', z_thet)
                    target_landmarks_rotated = self.rotateShape(target_landmarks)
                    temp = self.MSE(dict_landmarks, target_landmarks_rotated) 
                    if temp < mseMin:
                        mseMin = temp
                        best = [x_thet, y_thet, z_thet]
        return mseMin, np.array(best)
    def prune(self, A, B):
        # is dict_pose, B is new pose
        n_ldmks_A = A.shape[0]
        n_ldmks_B = B.shape[0]
        mseList = []
        thetaList = []
        if n_ldmks_A < n_ldmks_B:
            # prune B -- gives subset of B in shape of A
            combos = list(combinations(np.arange(n_ldmks_B), n_ldmks_A))
            for c in combos:
                subset = B[c, :]
                mse, theta = self.search(A, subset)
                mseList.append(mse)
                thetaList.append(theta)
            best_result_idx = np.argmin(mseList)
            return mseList[best_result_idx], np.array(thetaList)[best_result_idx, :]
        if n_ldmks_B < n_ldmks_A:
            # prune A -- gives subset of A in shape of B
            combos = list(combinations(np.arange(n_ldmks_A), n_ldmks_B))
            for c in combos:
                subset = A[c, :]
                mse, theta = self.search(subset, B)
                mseList.append(mse)
                thetaList.append(theta)

            best_result_idx = np.argmin(mseList)
            return mseList[best_result_idx], np.array(thetaList)[best_result_idx, :]
            



landmarks_all = {}
landmarks_base = {}
for file in os.listdir(os.path.join(os.getcwd(), 'landmarks_from_smpy')):
    if file != '.DS_Store':
        f_num = file[:file.index('.')]
        landmark_vals_coords = np.genfromtxt(os.path.join(os.getcwd(), 'landmarks_from_smpy', file), delimiter=',')
        parts = np.hsplit(landmark_vals_coords, [4])
        coords = parts[1]
        if f_num == '4' or f_num == '13' or f_num=='17' or f_num=='99':
            landmarks_base[str(f_num)] = coords
       
        landmarks_all[str(f_num)] = coords

comps_path2 = os.path.join(os.getcwd(), 'compsv2')
if not os.path.exists(comps_path2):
    os.mkdir(comps_path2)

for k in landmarks_all.keys():
    pc = ProcrustesComp()
    print(f'Finding match for {k}')
    
    x = pc.findBestMatch(landmarks_base, landmarks_all[k])
    array2 = np.zeros((len(x),4))
    
    for i in range(len(x)):
        print(f'x{i} {x[i]}')
        array2[i] = np.concatenate((np.array(x[i][0]).reshape(1,1), np.array([x[i][1][j] for j in range(3)]).reshape(3,1))).reshape(4,)
    
    dest = str(k) + ".csv" 
    np.savetxt(os.path.join(comps_path2, dest), array2, delimiter=',')

for k in landmarks_all.keys():
    dest = str(k) + ".csv"
    comps = np.genfromtxt(os.path.join(comps_path2, dest), delimiter=',')

    best = np.argmin(comps[:, 0])
    pose = list(landmarks_base.keys())[best]
    print(f'The closest match for file {k} was pose {pose}')