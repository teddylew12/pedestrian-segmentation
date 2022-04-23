import numpy as np
from LaPy.lapy.TriaIO import import_off
from LaPy.lapy import Solver
from scipy.io import savemat
from pathlib import Path
from copy import deepcopy

class Segmenter:
    def __init__(self, off_filename):
        self.mesh = import_off(off_filename)
        self.vertices = self.mesh.v
        self.triangles = self.mesh.t
        self.num_verts = self.vertices.shape[0]
        self.num_tris = self.triangles.shape[0]
        self.eigenvalues, self.eigenfunctions = Solver(self.mesh).eigs()
        self.efunc_one = self.eigenfunctions[:, 1]
        self.all_mins = []
        self.all_maxes = []
        self.all_saddles = []
        self.all_extremum = []

        self.selected_mins = []
        self.selected_maxes = []
        self.selected_saddles = []
        self.selected_extremum = []

        self.persistance_pairs = []
        self.noise_removal_steps = {}
        self.colors = np.zeros((self.num_tris, 1))
        self.allcolors = np.zeros((self.num_tris,1))
        self.spectral_colors = {}

        ### Visualization Tools
        self.animate_colors = False
        assert self.mesh.is_oriented()

    def run_all(self,name):
        self.get_min_max_saddles()
        self.construct_MS_complex()
        self.topological_noise_removal()
        threshold = .05
        step = max(self.noise_removal_steps.keys())-3
        self.segment(step, threshold)
        if name:
            self.save_to_matlab(name)
        landmark_dict = self.get_landmarks()
        print(landmark_dict)

    def get_min_max_saddles(self):
        for vtx in range(self.num_verts):
            # For each vertex, get neighbors and see if its a min,max or saddle
            adj_idx_efunc_vals = self.get_adjacent_vtx_efunc_vals(vtx)

            if self.efunc_one[vtx] >= np.max(adj_idx_efunc_vals):
                self.all_maxes.append(vtx)
            elif self.efunc_one[vtx] <= np.min(adj_idx_efunc_vals):
                self.all_mins.append(vtx)
            elif self.is_saddle(self.efunc_one[vtx], adj_idx_efunc_vals):
                self.all_saddles.append(vtx)
            else:
                pass
        self.all_extremum = self.all_maxes + self.all_mins

    def construct_MS_complex(self):
        pp = []
        for saddle in self.all_saddles:
            for ccw in self.get_ccw_neighbor_vtxs(saddle):
                if self.efunc_one[ccw] >= self.efunc_one[saddle]:
                    path = self.find_path(ccw, "up")
                    path.insert(0, saddle)
                    pp.append((saddle, path[-1], np.abs(self.efunc_one[saddle] - self.efunc_one[path[-1]])))

                else:
                    path = self.find_path(ccw, "down")
                    path.insert(0, saddle)
                    pp.append((saddle, path[-1], np.abs(self.efunc_one[saddle] - self.efunc_one[path[-1]])))
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

    def segment(self, optimal_noise_removal_step, level_set_threshold):
        '''
        Current hyperparameters are
        optimal_noise_removal =10
        level_set_threshold = .05
        but hopefully can be made automatic
        :param optimal_noise_removal_step:
        :param level_set_threshold:
        :return:
        '''
        persistance_pairs = self.noise_removal_steps[optimal_noise_removal_step]
        self.selected_saddles = [int(val) for val in np.unique(persistance_pairs[:, 0])]
        self.selected_mins = [int(val) for val in np.unique(persistance_pairs[:, 1]) if val in self.all_maxes]
        self.selected_maxes = [int(val) for val in np.unique(persistance_pairs[:, 1]) if val in self.all_mins]
        self.selected_extremum = self.selected_mins + self.selected_maxes
        boundary_triangles = self.get_boundary_triangles(level_set_threshold)
        self.color_to_boundary(boundary_triangles)

    def get_boundary_triangles(self, level_set_threshold):
        all_boundary_tris = []
        for saddle in self.selected_saddles:
            level_set_tris = self.findLevelSetOnMesh(saddle, level_set_threshold)
            all_boundary_tris.extend(list(level_set_tris))
        return list(np.unique(np.array(all_boundary_tris)))

    def findLevelSetOnMesh(self, saddle, epsilon):
        target = self.efunc_one[saddle]
        mask = np.logical_and(self.efunc_one > target - epsilon, self.efunc_one < target + epsilon)
        vtxs = np.where(mask == True)[0]

        return self.get_level_set_tris(vtxs)

    def get_level_set_tris(self, level_set_vtxs):
        all_level_set_tris = []
        # Vtx here is a 1x3 array of vals
        for vtx in level_set_vtxs:
            # Get row of vertices index
            all_level_set_tris.extend(list(np.unique(np.where(self.triangles == vtx)[0])))
        return np.unique(np.array(all_level_set_tris))

    def pick_correct_noise_removal_step(self, target_segmentations):
        '''
        Not implemented yet...
        Take self.noise_removal_steps...
        ASSUMES BOUNDARIES HAVE CLOSED LOOP
        work backwards through the keys, run segmentation and see how many unique colors you get
        Stop when you reach the correct amount (or you go above it if it skips)
        :param target_segmentations:
        :return:
        '''

    def get_optimal_level_set_threshold(self):
        '''
        Not implemented yet
        Start with 0.00 and add .05 until level sets actually form closed loops
        ASSUMES USING CORRECT NOISE REMOVAL STEP
        Run flooding algo with # number of extremum,
        check it we get that # of unique colors + 2 (unsegmented_middle and boundary)
        :return:
        '''

    def get_adjacent_vtx_efunc_vals(self, vtx):
        ccw_neighbors = self.get_ccw_neighbor_vtxs(vtx)
        return self.efunc_one[ccw_neighbors]

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
            adj_idx_efunc_vals = self.efunc_one[ccw_neighbors]
            # Get next vertex as the one with highest efunc value
            next_vtx = ccw_neighbors[compare(adj_idx_efunc_vals)]
            path.append(next_vtx)
            current_vtx = next_vtx
        return path

    def get_paths(self, saddle, extremum):
        '''
        Not implemented right now bc not needed except for visuals but basic idea is
        Figure out if extremum is min/max
        Get all neigbors that are lower/higher that starting point
        Run find path on them
        Save shortest of the paths
        :param saddle:
        :param extremum:
        :return:
        '''
        pass

    def color_to_boundary(self, boundary_tris):
        '''
        Breadth first search, marking seen vertices
        Assumes boundary tris segment mesh fully, otherwise the whole mesh will be explored in one call

        :param starting_vtx:
        :return:
        '''
        # Color Tracker
        # Color all boundaries
        step_ct = 0
        boundary_color = len(self.selected_extremum) + 1
        self.colors[boundary_tris] = 1
        self.spectral_colors[(np.inf, np.inf, np.inf, np.inf)] = 0
        if self.animate_colors:
            self.allcolors = deepcopy(self.colors)
        for i, crit in enumerate(self.selected_extremum):
            self.spectral_colors[self.get_spectral_embedding(crit, 4)] = (i + 1)/boundary_color
            color = (i + 1) /boundary_color
            queue = []
            visited = []
            # Get triangles that include starting vtx
            queue.extend(list(np.where(self.triangles == crit)[0]))
            # Loop until all possible triangles have been evaluated
            while queue:
                # Get unvisited triangle
                tri = queue.pop(0)
                if tri not in visited and tri not in boundary_tris:
                    # Visit triangle
                    visited.append(tri)
                    # Color triangle
                    self.colors[tri] = color
                    # Get neighboring triangles
                    for candidate_tri in self.get_neighbor_tris(tri):
                        # Add to queue if unvisited and not a boundary triangle
                        if candidate_tri not in visited and candidate_tri not in boundary_tris:
                            queue.append(candidate_tri)
                step_ct+=1
                if self.animate_colors:
                    if (step_ct%100)==0:
                        self.allcolors = np.concatenate((self.allcolors,self.colors),axis=1)


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

    def get_spectral_embedding(self, vtx, num_efuncs):
        '''
        Get spectral embedding of first num_efuncs eigenfunctions
        :param vtx:
        :param num_efuncs:
        :return:
        '''
        evals = self.eigenvalues[1:num_efuncs + 1]
        efuncs = self.eigenfunctions[vtx, 1:num_efuncs + 1]
        return tuple([efuncs[i] / evals[i] for i in range(num_efuncs)])

    def get_landmarks(self):
        '''
        Needs some work
        Outputs a consistent landmark matrix
        :return:
        '''
        crit_means = {}
        for spectral_embedding, color in self.spectral_colors.items():
            tris = self.triangles[(self.colors == color).flatten()]
            crit_means[spectral_embedding] = np.mean(self.vertices[tris.flatten()], axis=0)
        return crit_means

    def save_to_matlab(self, name):
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
        output_data["eigenfunction"] = self.efunc_one
        output_data["eigenfunctions"]=self.eigenfunctions

        output_data["colors"] = self.colors
        if self.animate_colors:
            output_data["allcolors"]=self.allcolors
        output_file = Path.cwd() / (name + ".mat")
        savemat(output_file, output_data)
    def display_mesh(self,save_name):
        from LaPy.lapy.Plot import plot_tria_mesh
        plot_tria_mesh(self.mesh,export_png=f"{save_name}.png")
model = "./off/14.off"
segmenter = Segmenter(model)
segmenter.animate_colors = True
segmenter.run_all("allcolors2")
