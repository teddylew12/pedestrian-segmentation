import numpy as np
from LaPy.lapy.TriaIO import import_off
from LaPy.lapy import Solver
from scipy.io import savemat
from pathlib import Path
from copy import deepcopy
from sklearn.preprocessing import MinMaxScaler
class CriticalPoint:
    def __init__(self, vertex,embedding,color = -1):
        self.vertex =vertex
        self.color = color
        self.embedding = embedding

class Segmenter:
    def __init__(self, off_filename):
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
        self.colors = np.zeros((self.num_tris, 1))
        self.allcolors = np.zeros((self.num_tris,1))
        self.spectral_colors = {}
        self.crit_embeddings = {}
        self.crits=[]

        ### Visualization Tools
        self.isbasemesh = True
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
        #print(landmark_dict)
        return landmark_dict
    def flip_eigenfunctions(self, base_mesh):
        '''
        Pseudocode
        For Each eigenfunction of the current mesh
            Find max location of that efunc and of the flipped efunc... find max location of that efunc on the other mesh
            if the vtx distance is closer with the other efunc and within a certain distance (if they get switched somehow)
            then flip it!

        :param base_mesh:
        :return:
        '''
        for i in range(self.eigenfunctions.shape[1]):
            if i==0:
                continue
            max_on_base_shape = base_mesh.vertices[np.argmax(base_mesh.eigenfunctions[:,i])]
            max_on_current_shape = self.vertices[np.argmax(self.eigenfunctions[:,i])]
            flipped_max_on_current_shape = self.vertices[np.argmax(-1*self.eigenfunctions[:,i])]

            main_dist = np.linalg.norm(max_on_base_shape-max_on_current_shape)
            flipped_dist = np.linalg.norm(max_on_base_shape-flipped_max_on_current_shape)
            if flipped_dist<main_dist:
                self.eigenfunctions[:,i] = -1 * self.eigenfunctions[:,i]
        self.eigenfunctions[:,2] = -1 * self.eigenfunctions[:,2]



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
        self.boundary_triangles = self.get_boundary_triangles(level_set_threshold)
        if self.isbasemesh:
            self.color_to_boundary()


    def get_boundary_triangles(self, level_set_threshold):
        all_boundary_tris = []
        for saddle in self.selected_saddles:
            level_set_tris = self.findLevelSetOnMesh(saddle, level_set_threshold)
            all_boundary_tris.extend(list(level_set_tris))
        return list(np.unique(np.array(all_boundary_tris)))

    def findLevelSetOnMesh(self, saddle, epsilon):
        target = self.eigenfunctions[:,1][saddle]
        mask = np.logical_and(self.eigenfunctions[:,1] > target - epsilon, self.eigenfunctions[:,1] < target + epsilon)
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

    def color_to_boundary(self, ):
        '''
        Breadth first search, marking seen vertices
        Assumes boundary tris segment mesh fully, otherwise the whole mesh will be explored in one call

        :param starting_vtx:
        :return:
        '''
        # Color Tracker
        # Color all boundaries
        num_embedding_dims = 4
        step_ct = 0
        boundary_color = len(self.selected_extremum) + 1
        self.colors[self.boundary_triangles] = 1
        self.spectral_colors[tuple([100000]*num_embedding_dims)] = 0
        if self.animate_colors:
            self.allcolors = deepcopy(self.colors)
        for i, crit in enumerate(self.selected_extremum):
            embedding = self.get_spectral_embedding(crit, num_embedding_dims)
            self.crit_embeddings[embedding] = crit
            self.spectral_colors[embedding] = (i + 1)/boundary_color
            color = (i + 1) /boundary_color
            queue = []
            visited = []
            # Get triangles that include starting vtx
            queue.extend(list(np.where(self.triangles == crit)[0]))
            # Loop until all possible triangles have been evaluated
            while queue:
                # Get unvisited triangle
                tri = queue.pop(0)
                if tri not in visited and tri not in self.boundary_triangles:
                    # Visit triangle
                    visited.append(tri)
                    # Color triangle
                    self.colors[tri] = color
                    # Get neighboring triangles
                    for candidate_tri in self.get_neighbor_tris(tri):
                        # Add to queue if unvisited and not a boundary triangle
                        if candidate_tri not in visited and candidate_tri not in self.boundary_triangles:
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
        output_data["eigenfunction"] = self.eigenfunctions[:,1]
        output_data["eigenfunctions"]=self.eigenfunctions

        output_data["colors"] = self.colors
        if self.animate_colors:
            output_data["allcolors"]=self.allcolors
        output_file = Path.cwd()/"runs" / (name + ".mat")
        savemat(output_file, output_data)
    def display_mesh(self,save_name):
        from LaPy.lapy.Plot import plot_tria_mesh
        plot_tria_mesh(self.mesh,export_png=f"{save_name}.png")

    def get_critical_points(self):
        for crit in self.crit_embeddings:
            try:
                self.crits.append(CriticalPoint(self.crit_embeddings[crit],crit,color =self.spectral_colors[crit]))
            except:
                self.crits.append(CriticalPoint(self.crit_embeddings[crit], crit))
    def get_crit_embeddings(self):
        for crit in self.selected_extremum:
            embedding = self.get_spectral_embedding(crit,4)
            self.crit_embeddings[embedding] = crit

    def color_to_boundary_matching(self, color_dict,fname):
        '''
        Breadth first search, marking seen vertices
        Assumes boundary tris segment mesh fully, otherwise the whole mesh will be explored in one call

        :param starting_vtx:
        :return:
        '''
        # Color Tracker
        # Color all boundaries
        self.colors = np.zeros((self.num_tris,1))
        step_ct = 0
        self.colors[self.boundary_triangles] = 1
        if self.animate_colors:
            self.allcolors = deepcopy(self.colors)
        for crit, color in color_dict.items():
            queue = []
            visited = []
            # Get triangles that include starting vtx
            queue.extend(list(np.where(self.triangles == crit)[0]))
            # Loop until all possible triangles have been evaluated
            while queue:
                # Get unvisited triangle
                tri = queue.pop(0)
                if tri not in visited and tri not in self.boundary_triangles:
                    # Visit triangle
                    visited.append(tri)
                    # Color triangle
                    self.colors[tri] = color
                    # Get neighboring triangles
                    for candidate_tri in self.get_neighbor_tris(tri):
                        # Add to queue if unvisited and not a boundary triangle
                        if candidate_tri not in visited and candidate_tri not in self.boundary_triangles:
                            queue.append(candidate_tri)
                step_ct+=1
                if self.animate_colors:
                    if (step_ct%100)==0:
                        self.allcolors = np.concatenate((self.allcolors,self.colors),axis=1)
        self.save_to_matlab(fname)
    def match_segment_colors(self,base_mesh_colors, fname):
        self.get_crit_embeddings()
        self.get_critical_points()
        color_correspondances = {}
        base_embeddings_array = np.array([list(embedding)for embedding in base_mesh_colors.keys()])
        for crit in self.crits:
            candidate = base_embeddings_array - np.array(list(crit.embedding)).reshape(1, 4)
            norms = np.linalg.norm(candidate, axis=1)
            closest = np.argmin(norms)
            closest_embedding = base_embeddings_array[closest]
            color_correspondances[crit.vertex] = base_mesh_colors[tuple(closest_embedding)]
            self.color_to_boundary_matching(color_correspondances,fname)
model4 = "./off/4.off"
base_mesh = Segmenter(model4)
lmdict1=base_mesh.run_all("base_mesh")
base_mesh.get_critical_points()


model13 = "./off/13.off"
segmenter13 = Segmenter(model13)
segmenter13.isbasemesh=False
segmenter13.flip_eigenfunctions(base_mesh)
lmdict2= segmenter13.run_all("")
segmenter13.match_segment_colors(base_mesh.spectral_colors, "matching")
