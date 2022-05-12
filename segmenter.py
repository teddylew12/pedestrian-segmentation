from LaPy.lapy.TriaIO import import_off
from LaPy.lapy import Solver
import numpy as np
from scipy.io import savemat
from sklearn.preprocessing import MinMaxScaler
from proctrustes import ProcrustesComp
import os

class Segmenter:
    def __init__(self, off_filename, results_dir):
        self.results_dir = results_dir
        print(off_filename)
        self.mesh = import_off(off_filename)
        self.vertices = self.mesh.v
        self.triangles = self.mesh.t
        self.num_verts = self.vertices.shape[0]
        self.num_tris = self.triangles.shape[0]
        self.eigenvalues, self.eigenfunctions = Solver(self.mesh).eigs()
        scaler = MinMaxScaler(feature_range=(-1,1))
        self.eigenfunctions = scaler.fit_transform(self.eigenfunctions)

        self.all_mins = []
        self.all_maxes = []
        self.all_saddles = []
        self.all_extremum = []

        self.selected_mins = []
        self.selected_maxes = []
        self.selected_saddles = []
        self.selected_extremum = []
        self.boundary_triangles = []

        self.persistance_pairs = []
        self.noise_removal_steps = {}
        self.full_landmark_matrix = []
        self.landmark_points_matrix = []

        self.level_set_threshold = .05
        self.colors = np.zeros((self.num_tris, 1))

        assert self.mesh.is_oriented()

    def save_initial(self):
        name = "initial_eigenfunctions.mat"
        output_data = {"eigenfunctions": self.eigenfunctions, "triangles": self.triangles, "vertices":self.vertices}
        output_file = self.results_dir / name
        savemat(output_file, output_data)

    def flip_eigenfunctions(self):
        self.save_initial()
        print("Open matlab and run compare_eigenfunctions to visualize base and current eigenfunctions!")
        print("You will need to update the filepath for the non-base mesh")
        for i in range(1,5):
            inp = input(f"Flip eigenfunction {i}?\nType y for Yes, n for No:")
            if inp == "y":
                print(f"Flipping eigenfunction {i}")
                self.eigenfunctions[:,i] *= -1



    def get_min_max_saddles(self):
        for vtx in range(self.num_verts):
            # For each vertex, get neighbors and see if its a min,max or saddle
            adj_idx_efunc_vals = self.get_adjacent_vtx_efunc_vals(vtx)

            if self.eigenfunctions[:,1][vtx] >= np.max(adj_idx_efunc_vals):
                self.all_maxes.append(vtx)
            elif self.eigenfunctions[:,1][vtx] <= np.min(adj_idx_efunc_vals):
                self.all_mins.append(vtx)
            elif self.is_saddle(self.eigenfunctions[:,1][vtx], adj_idx_efunc_vals):
                self.all_saddles.append(vtx)
            else:
                pass
        self.all_extremum = self.all_maxes + self.all_mins

    def construct_MS_complex(self):
        pp = []
        for saddle in self.all_saddles:
            for ccw in self.get_ccw_neighbor_vtxs(saddle):
                if self.eigenfunctions[:,1][ccw] >= self.eigenfunctions[:,1][saddle]:
                    path = self.find_path(ccw, "up")
                    path.insert(0, saddle)
                    pp.append((saddle, path[-1], np.abs(self.eigenfunctions[:,1][saddle] - self.eigenfunctions[:,1][path[-1]])))

                else:
                    path = self.find_path(ccw, "down")
                    path.insert(0, saddle)
                    pp.append((saddle, path[-1], np.abs(self.eigenfunctions[:,1][saddle] - self.eigenfunctions[:,1][path[-1]])))
        pp = np.array([[a, b, c] for a, b, c in pp])
        pp = np.unique(pp, axis=0)
        self.persistance_pairs = pp

    def topological_noise_removal(self):
        self.persistance_pairs = self.persistance_pairs[np.argsort(self.persistance_pairs[:, 2])]
        step_ct = 0
        self.noise_removal_steps[step_ct] = self.persistance_pairs
        while self.persistance_pairs.size > 0:
            self.persistance_pairs = self.persistance_pairs[np.argsort(self.persistance_pairs[:, 2])]
            step_ct += 1
            saddle, extremum, _ = self.persistance_pairs[0]
            self.persistance_pairs = np.delete(self.persistance_pairs, 0, 0)
            other_pts_with_same_saddle = self.persistance_pairs[self.persistance_pairs[:, 0] == saddle]
            if extremum in self.all_mins:
                for ext in other_pts_with_same_saddle[:, 1]:
                    if ext in self.all_mins:
                        # Anywhere where extremum is should be replaced with ext
                        self.persistance_pairs[self.persistance_pairs == extremum] = ext
                        self.persistance_pairs = self.persistance_pairs[self.persistance_pairs[:, 0] != saddle]
                        break
            else:
                for ext in other_pts_with_same_saddle[:, 1]:
                    if ext in self.all_maxes:
                        # Anywhere where extremum is should be replaced with ext
                        self.persistance_pairs[self.persistance_pairs == extremum] = ext
                        self.persistance_pairs = self.persistance_pairs[self.persistance_pairs[:, 0] != saddle]
                        break
            self.noise_removal_steps[step_ct] = self.persistance_pairs

    def get_selected_extremum(self, step):
        persistance_pairs = self.noise_removal_steps[step]
        self.selected_saddles = [int(val) for val in np.unique(persistance_pairs[:, 0])]
        self.selected_mins = [int(val) for val in np.unique(persistance_pairs[:, 1]) if val in self.all_maxes]
        self.selected_maxes = [int(val) for val in np.unique(persistance_pairs[:, 1]) if val in self.all_mins]
        self.selected_extremum = self.selected_mins + self.selected_maxes

    def get_landmark_matrix(self):
        eps = .075
        landmark = np.zeros((len(self.selected_extremum),8))
        landmark[:,:4] = self.eigenfunctions[self.selected_extremum,1:5]
        landmark[:,4:7]=self.vertices[self.selected_extremum]
        landmark[:,7]=self.selected_extremum
        landmark = landmark[np.argsort(landmark[:,0])]
        for row in range(landmark.shape[0]-1):
            col = 0
            if np.abs(landmark[row,col]-landmark[row+1,col])<eps and col <=4:
                while np.abs(landmark[row,col]-landmark[row+1,col])<eps:
                    col+=1
                if landmark[row,col] > landmark[row+1,col]:
                    landmark[[row,row+1]] = landmark[[row+1,row]]
        self.full_landmark_matrix = landmark
        self.landmark_points_matrix = landmark[:, 4:7]

    def get_boundary_triangles(self):
        all_boundary_tris = []
        for saddle in self.selected_saddles:
            level_set_tris = self.findLevelSetOnMesh(saddle)
            all_boundary_tris.extend(list(level_set_tris))
        return list(np.unique(np.array(all_boundary_tris)))

    def findLevelSetOnMesh(self, saddle):
        epsilon = self.level_set_threshold
        target = self.eigenfunctions[:,1][saddle]
        mask = np.logical_and(self.eigenfunctions[:,1] > target - epsilon, self.eigenfunctions[:,1] < target + epsilon)
        level_set_vtxs = np.where(mask == True)[0]

        all_level_set_tris = []
        # Vtx here is a 1x3 array of vals
        for vtx in level_set_vtxs:
            # Get row of vertices index
            all_level_set_tris.extend(list(np.unique(np.where(self.triangles == vtx)[0])))
        return np.unique(np.array(all_level_set_tris))


    def get_adjacent_vtx_efunc_vals(self, vtx):
        ccw_neighbors = self.get_ccw_neighbor_vtxs(vtx)
        return self.eigenfunctions[:,1][ccw_neighbors]

    def get_ccw_neighbor_vtxs(self, vtx):
        # Mesh is oriented, so each triangle is listed in ccw order
        # ALGO
        # 1. Get all neighboring triangles
        # 2. Pick one triangle (doesnt matter which one), get first 2 neighbors in ccw by cycling through
        # i.e. if the target vertex is 25 and the target triangle is 26 25 14, we get 14 and then 26
        # 3. Find the next triangle, which is the triangle that contains the vertex and the last added neighbor
        # 4. Repeat around the whole neighborhood!
        neighbors = []
        # Get rows, cols where the target vertex exists in the triangle
        vtx_locations = np.argwhere(self.triangles == vtx)
        # Get triangles that contain the vertex
        adj_triangles = self.triangles[vtx_locations[:, 0]]
        # Get first triangle
        first_r, first_c = vtx_locations[0]
        # Get first 2 neighbors in ccw manner
        neighbors.append(self.triangles[first_r, (first_c + 1) % 3])
        neighbors.append(self.triangles[first_r, (first_c + 2) % 3])
        # Keep track of last added neighbor
        prev = self.triangles[first_r, (first_c + 2) % 3]
        # Remove triangle from list of possible tris
        adj_triangles = np.delete(adj_triangles, 0, 0)
        num_tris = adj_triangles.shape[0]
        # Repeat for other triangles
        for _ in range(num_tris - 1):
            for i, tri in enumerate(adj_triangles):
                # Find triangle that shares an edge with last processed triangle
                if prev in tri:
                    # Get new vertex
                    new_element = set(tri) - {prev, vtx}
                    # Add to neighbors
                    neighbors.append(list(new_element)[0])
                    # Update prev and available triangles
                    prev = list(new_element)[0]
                    adj_triangles = np.delete(adj_triangles, i, 0)
                    break
        return neighbors

    def is_saddle(self, center_vtx_val, adj_vals):
        '''
         A saddle is a vertex i where the value of hj − hi at the neighbors j
         switches sign at least 4 time
        :param val:
        :param adj_vals:
        :return:
        '''
        # Add first value to end to test all ccw pairs of verticles
        adj_vals = np.append(adj_vals, adj_vals[0])
        # Get val of hj-hi
        func_diffs = adj_vals - center_vtx_val
        # Find where it crosses 0
        # @Oscar this might be buggy idk
        # Professor said you can also calculate gaussian curvature and see if its negative
        # So maybe we also try to implement that
        zero_crossings = np.where(np.diff(np.sign(func_diffs)))[0]
        return len(zero_crossings) >= 4

    def find_path(self, saddle_vtx, direction):
        '''
        At each vertex, compute gradient and pick edge that is more in direction of gradient
        then go to that
        :param saddle_vtx:
        :return:
        '''
        path = [saddle_vtx]
        current_vtx = saddle_vtx
        # Keep going until we end at a max vertex
        if direction == "up":
            target = self.all_maxes
            compare = np.argmax
        else:
            target = self.all_mins
            compare = np.argmin
        while current_vtx not in target:
            # Get adjacent vertices
            ccw_neighbors = [neighbor for neighbor in self.get_ccw_neighbor_vtxs(current_vtx) if neighbor not in path]
            adj_idx_efunc_vals = self.eigenfunctions[:,1][ccw_neighbors]
            # Get next vertex as the one with highest efunc value
            next_vtx = ccw_neighbors[compare(adj_idx_efunc_vals)]
            path.append(next_vtx)
            current_vtx = next_vtx
        return path


    def color_mesh(self):
        '''
        Breadth first search, marking seen vertices
        Assumes boundary tris segment mesh fully, otherwise the whole mesh will be explored in one call
        :return:
        '''
        boundary_triangles = self.get_boundary_triangles()
        boundary_color = self.full_landmark_matrix.shape[0] + 1
        self.colors[boundary_triangles]=boundary_color
        for i, row in enumerate(self.full_landmark_matrix):
            color = i + 1
            crit = row[-1]
            queue = []
            visited = []
            # Get triangles that include starting vtx
            queue.extend(list(np.where(self.triangles == crit)[0]))
            # Loop until all possible triangles have been evaluated
            while queue:
                # Get unvisited triangle
                tri = queue.pop(0)
                if tri not in visited and tri not in boundary_triangles:
                    # Visit triangle
                    visited.append(tri)
                    # Color triangle
                    self.colors[tri] = color
                    # Get neighboring triangles
                    for candidate_tri in self.get_neighbor_tris(tri):
                        # Add to queue if unvisited and not a boundary triangle
                        if candidate_tri not in visited and candidate_tri not in boundary_triangles:
                            queue.append(candidate_tri)

    def get_neighbor_tris(self, tri):
        all_tris = []
        # Get vertices in the target triangle
        for vtx in self.triangles[tri]:
            # Get all triangles that share the vertex
            neigbor_tris = np.where(self.triangles == vtx)[0]
            # Add to list of neighboring triangles
            all_tris.extend(list(neigbor_tris))
        # Remove duplicates
        return list(np.unique(np.array(all_tris)))

    def save_to_matlab(self):
        '''
        Save all calculated data for matlab viewing
        :param name:
        :return:
        '''

        output_data = {}
        output_data["all_maxes"] = self.all_maxes
        output_data["all_mins"] = self.all_mins
        output_data["all_saddles"] = self.all_saddles

        output_data["select_maxes"] = self.selected_maxes
        output_data["select_mins"] = self.selected_mins
        output_data["select_saddles"] = self.selected_saddles

        output_data["vertices"] = self.vertices
        output_data["triangles"] = self.triangles
        output_data["eigenfunctions"]=self.eigenfunctions

        output_data["colors"] = self.colors
        output_data["landmark_matrix"] = self.full_landmark_matrix

        output_file = self.results_dir / "base_mesh.mat"
        savemat(output_file, output_data)
    

    def get_closest_matching_pose(self):
        target_pose = self.landmark_points_matrix
       
        poses_dict = {}
        for file in os.listdir(os.path.join(os.getcwd(), 'pose_dictionary_landmarks')):
            if file != '.DS_Store':
                f_num = file[:file.index('.')]
                landmark_vals_coords = np.genfromtxt(os.path.join(os.getcwd(), 'pose_dictionary_landmarks', file), delimiter=',')
                parts = np.hsplit(landmark_vals_coords, [4])
                coords = parts[1]
                if f_num == '4' or f_num == '13' or f_num=='17' or f_num=='99':
                    poses_dict[str(f_num)] = coords
                else:
                    print('Unidentified landmark pose in file!: ' + str(f_num))
       
                poses_dict[str(f_num)] = coords
        
        pc = ProcrustesComp()
        
        matches = pc.findBestMatch(poses_dict, target_pose)
        mse = [matches[i][0] for i in range(len(matches))]
        best = np.argmin(mse)

        pose = list(poses_dict.keys())[best]
        if str(pose) == '4':
            return "running"
        elif str(pose) ==  '13':
            return "standing"
        elif str(pose) == '17':
            return "sitting"
        elif str(pose) == '99':
            return "laying down"
        else:
            print('hmmm matching didnt go right')
            return -1




