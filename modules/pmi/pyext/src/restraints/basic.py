"""@namespace IMP.pmi.restraints.basic
Some miscellaneous simple restraints.
"""

import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.pmi.tools
import IMP.pmi.restraints
import math


class ExternalBarrier(IMP.pmi.restraints.RestraintBase):
    """Keeps all structures inside a sphere."""

    def __init__(self, hierarchies, radius=10.0, resolution=10, weight=1.0,
                 center=None, label=None):
        """Setup external barrier restraint.
        @param hierarchies Can be one of the following inputs: IMP Hierarchy,
               PMI System/State/Molecule/TempResidue, or a list/set of them
        @param radius Size of external barrier
        @param resolution Select which resolutions to act upon
        @param weight Weight of restraint
        @param center Center of the external barrier
               (IMP.algebra.Vector3D object)
        @param label A unique label to be used in outputs and
                     particle/restraint names.
        """
        hiers = IMP.pmi.tools.input_adaptor(hierarchies, resolution,
                                            flatten=True)
        model = hiers[0].get_model()
        particles = [h.get_particle() for h in hiers]

        super().__init__(model, label=label, weight=weight)
        self.radius = radius

        if center is None:
            c3 = IMP.algebra.Vector3D(0, 0, 0)
        elif type(center) is IMP.algebra.Vector3D:
            c3 = center
        else:
            raise Exception(
                "%s: @param center must be an IMP.algebra.Vector3D object" % (
                    self.name))

        ub3 = IMP.core.HarmonicUpperBound(radius, 10.0)
        ss3 = IMP.core.DistanceToSingletonScore(ub3, c3)
        lsc = IMP.container.ListSingletonContainer(self.model)

        lsc.add(particles)
        r3 = IMP.container.SingletonsRestraint(ss3, lsc)
        self.rs.add_restraint(r3)


class DistanceRestraint(IMP.pmi.restraints.RestraintBase):
    """A simple distance restraint"""

    def __init__(self, root_hier, tuple_selection1, tuple_selection2,
                 distancemin=0, distancemax=100, resolution=1.0, kappa=1.0,
                 label=None, weight=1.):
        """Setup distance restraint.
        @param root_hier The hierarchy to select from
        @param tuple_selection1 (resnum, resnum, molecule name, copy
               number (=0))
        @param tuple_selection2 (resnum, resnum, molecule name, copy
               number (=0))
        @param distancemin The minimum dist
        @param distancemax The maximum dist
        @param resolution For selecting particles
        @param kappa The harmonic parameter
        @param label A unique label to be used in outputs and
                     particle/restraint names
        @param weight Weight of restraint
        @note Pass the same resnum twice to each tuple_selection. Optionally
              add a copy number.
        """
        ts1 = IMP.core.HarmonicUpperBound(distancemax, kappa)
        ts2 = IMP.core.HarmonicLowerBound(distancemin, kappa)

        model = root_hier.get_model()
        copy_num1 = 0
        if len(tuple_selection1) > 3:
            copy_num1 = tuple_selection1[3]
        copy_num2 = 0
        if len(tuple_selection2) > 3:
            copy_num2 = tuple_selection2[3]

        sel1 = IMP.atom.Selection(root_hier,
                                  resolution=resolution,
                                  molecule=tuple_selection1[2],
                                  residue_index=tuple_selection1[0],
                                  copy_index=copy_num1)
        particles1 = sel1.get_selected_particles()
        sel2 = IMP.atom.Selection(root_hier,
                                  resolution=resolution,
                                  molecule=tuple_selection2[2],
                                  residue_index=tuple_selection2[0],
                                  copy_index=copy_num2)
        particles2 = sel2.get_selected_particles()

        super().__init__(model, label=label, weight=weight)
        print(self.name)

        print("Created distance restraint between "
              "%s and %s" % (particles1[0].get_name(),
                             particles2[0].get_name()))

        if len(particles1) > 1 or len(particles2) > 1:
            raise ValueError("more than one particle selected")

        self.rs.add_restraint(
            IMP.core.DistanceRestraint(self.model, ts1,
                                       particles1[0],
                                       particles2[0]))
        self.rs.add_restraint(
            IMP.core.DistanceRestraint(self.model, ts2,
                                       particles1[0],
                                       particles2[0]))


class CylinderRestraint(IMP.Restraint):
    """Restrain particles within (or outside) a cylinder.
       The cylinder is aligned along the z-axis and with center x=y=0.
       Optionally, one can restrain the cylindrical angle
    """
    import math

    def __init__(self, m, objects, resolution, radius, mintheta=None,
                 maxtheta=None, repulsive=False, label='None'):
        '''
        @param objects PMI2 objects to restrain
        @param resolution the resolution you want the restraint to be applied
        @param radius the radius of the cylinder
        @param mintheta minimum cylindrical angle in degrees
        @param maxtheta maximum cylindrical angle in degrees
        @param repulsive If True, restrain the particles to be outside
               of the cylinder instead of inside
        @param label A unique label to be used in outputs and
               particle/restraint names
        '''
        IMP.Restraint.__init__(self, m, "CylinderRestraint %1%")
        self.radius = radius
        self.softness = 3.0
        self.softness_angle = 0.5
        self.plateau = 1e-10
        self.weight = 1.0
        self.m = m
        self.mintheta = mintheta
        self.maxtheta = maxtheta
        self.repulsive = repulsive
        hierarchies = IMP.pmi.tools.input_adaptor(objects,
                                                  resolution,
                                                  flatten=True)
        self.particles = [h.get_particle() for h in hierarchies]
        self.label = label

    def get_probability(self, p):
        xyz = IMP.core.XYZ(p)
        r = self.math.sqrt(xyz.get_x()**2+xyz.get_y()**2)
        argvalue = (r-self.radius) / self.softness
        if self.repulsive:
            argvalue = -argvalue
        prob = (1.0 - self.plateau) / (1.0 + self.math.exp(-argvalue))
        return prob

    def get_angle_probability(self, p):
        xyz = IMP.core.XYZ(p)
        angle = self.math.atan2(xyz.get_y(), xyz.get_x())*180.0/self.math.pi
        anglediff = (angle - self.maxtheta + 180 + 360) % 360 - 180
        argvalue1 = anglediff / self.softness_angle
        anglediff = (angle - self.mintheta + 180 + 360) % 360 - 180
        argvalue2 = -anglediff / self.softness_angle
        prob = ((1.0-self.plateau)
                / (1.0 + self.math.exp(-max(argvalue1, argvalue2))))
        return prob

    def unprotected_evaluate(self, da):
        s = 0.0
        for p in self.particles:
            s += -self.math.log(1.0-self.get_probability(p))
            if self.mintheta is not None and self.maxtheta is not None:
                s += -self.math.log(1.0-self.get_angle_probability(p))
        return s

    def do_get_inputs(self):
        return self.particles

    def add_to_model(self):
        IMP.pmi.tools.add_restraint_to_model(self.m, self)

    def get_output(self):
        output = {}
        score = self.weight * self.unprotected_evaluate(None)
        output["_TotalScore"] = str(score)
        output["CylinderRestraint_" + self.label] = str(score)
        return output


class BiStableDistanceRestraint(IMP.Restraint):
    '''Distance restraint with bistable potential
    Authors: G. Bouvier, R. Pellarin. Pasteur Institute.
    '''
    import numpy as np
    import math

    def __init__(self, m, p1, p2, dist1, dist2, sigma1, sigma2, weight1,
                 weight2):
        '''
        input two particles, the two equilibrium distances, their amplitudes,
        and their weights (populations)
        '''
        IMP.Restraint.__init__(self, m, "BiStableDistanceRestraint %1%")
        self.dist1 = dist1
        self.dist2 = dist2

        self.sigma1 = sigma1
        self.sigma2 = sigma2

        self.weight1 = weight1
        self.weight2 = weight2

        if self.weight1+self.weight2 != 1:
            raise ValueError("The sum of the weights must be one")

        self.d1 = IMP.core.XYZ(p1)
        self.d2 = IMP.core.XYZ(p2)
        self.particle_list = [p1, p2]

    def gaussian(self, x, mu, sig, w):
        return (w*self.np.exp(-self.np.power(x - mu, 2.)
                / (2 * self.np.power(sig, 2.))))

    def unprotected_evaluate(self, da):
        dist = IMP.core.get_distance(self.d1, self.d2)
        prob = self.gaussian(dist, self.dist1, self.sigma1, self.weight1) + \
            self.gaussian(dist, self.dist2, self.sigma2, self.weight2)
        return -self.math.log(prob)

    def do_get_inputs(self):
        return self.particle_list


class DistanceToPointRestraint(IMP.pmi.restraints.RestraintBase):
    """Anchor a particle to a specific coordinate."""

    def __init__(self, root_hier, tuple_selection,
                 anchor_point=IMP.algebra.Vector3D(0, 0, 0),
                 radius=10.0, kappa=10.0, resolution=1.0, weight=1.0,
                 label=None):
        """Setup distance restraint.
        @param root_hier The hierarchy to select from
        @param tuple_selection (resnum, resnum, molecule name,
               copy number (=0))
        @param anchor_point Point to which to restrain particle
               (IMP.algebra.Vector3D object)
        @param radius Maximum distance the particle can move from the
               fixed point before it is restrained
        @param kappa Strength of the harmonic restraint for point-particle
               distance greater than the radius
        @param resolution For selecting a particle
        @param weight Weight of restraint
        @param label A unique label to be used in outputs and
                     particle/restraint names
        @note Pass the same resnum twice to each tuple_selection. Optionally
              add a copy number
        """
        model = root_hier.get_model()
        copy_num1 = 0
        if len(tuple_selection) > 3:
            copy_num1 = tuple_selection[3]

        sel1 = IMP.atom.Selection(root_hier,
                                  resolution=resolution,
                                  molecule=tuple_selection[2],
                                  residue_index=tuple_selection[0],
                                  copy_index=copy_num1)
        ps = sel1.get_selected_particles()
        if len(ps) > 1:
            raise ValueError("More than one particle selected")

        super().__init__(model, label=label, weight=weight)
        self.radius = radius

        ub3 = IMP.core.HarmonicUpperBound(self.radius, kappa)
        if anchor_point is None:
            c3 = IMP.algebra.Vector3D(0, 0, 0)
        elif isinstance(anchor_point, IMP.algebra.Vector3D):
            c3 = anchor_point
        else:
            raise TypeError("anchor_point must be an algebra.Vector3D object")
        ss3 = IMP.core.DistanceToSingletonScore(ub3, c3)

        lsc = IMP.container.ListSingletonContainer(self.model)
        lsc.add(ps)

        r3 = IMP.container.SingletonsRestraint(ss3, lsc)
        self.rs.add_restraint(r3)

        print("\n%s: Created distance_to_point_restraint between "
              "%s and %s" % (self.name, ps[0].get_name(), c3))


class MembraneRestraint(IMP.pmi.restraints.RestraintBase):
    """Restrain particles to be above, below, or inside a planar membrane.
       The membrane is defined to lie on the xy plane with a given z
       coordinate and thickness, and particles are restrained (by their
       z coordinates) with a simple sigmoid score.
    """

    def __init__(self, hier, objects_above=None, objects_inside=None,
                 objects_below=None, center=0.0, thickness=30.0,
                 softness=3.0, plateau=0.0000000001, resolution=1,
                 weight=1.0, label=None):
        """Setup the restraint.

        @param objects_inside list or tuples of objects in membrane
               (e.g. ['p1', (10, 30,'p2')])
        @param objects_above list or tuples of objects above membrane
        @param objects_below list or tuples of objects below membrane
        @param thickness Thickness of the membrane along the z-axis
        @param softness Softness of the limiter in the sigmoid function
        @param plateau Parameter to set the probability (=1- plateau))
               at the plateau phase of the sigmoid
        @param weight Weight of restraint
        @param label A unique label to be used in outputs and
                     particle/restraint names.
        """

        self.hier = hier
        model = self.hier.get_model()

        super().__init__(
            model, name="MembraneRestraint", label=label, weight=weight)

        self.center = center
        self.thickness = thickness
        self.softness = softness
        self.plateau = plateau
        self.linear = 0.02
        self.resolution = resolution

        # Create nuisance particle
        p = IMP.Particle(model)
        z_center = IMP.isd.Nuisance.setup_particle(p)
        z_center.set_nuisance(self.center)

        # Setup restraint
        mr = IMP.pmi.MembraneRestraint(model,
                                       z_center.get_particle_index(),
                                       self.thickness,
                                       self.softness,
                                       self.plateau,
                                       self.linear)

        # Particles above
        if objects_above:
            for obj in objects_above:
                if isinstance(obj, tuple):
                    self.particles_above = self._select_from_tuple(obj)

                elif isinstance(obj, str):
                    self.particles_above = self._select_from_string(obj)
                mr.add_particles_above(self.particles_above)

        # Particles inside
        if objects_inside:
            for obj in objects_inside:
                if isinstance(obj, tuple):
                    self.particles_inside = self._select_from_tuple(obj)

                elif isinstance(obj, str):
                    self.particles_inside = self._select_from_string(obj)
                mr.add_particles_inside(self.particles_inside)

        # Particles below
        if objects_below:
            for obj in objects_below:
                if isinstance(obj, tuple):
                    self.particles_below = self._select_from_tuple(obj)

                elif isinstance(obj, str):
                    self.particles_below = self._select_from_string(obj)
                mr.add_particles_below(self.particles_below)

        self.rs.add_restraint(mr)

    def get_particles_above(self):
        return self.particles_above

    def get_particles_inside(self):
        return self.particles_inside

    def get_particles_below(self):
        return self.particles_below

    def _select_from_tuple(self, obj):
        particles = IMP.atom.Selection(
            self.hier, molecule=obj[2],
            residue_indexes=range(obj[0], obj[1]+1, 1),
            resolution=self.resolution).get_selected_particles()

        return particles

    def _select_from_string(self, obj):
        particles = IMP.atom.Selection(
            self.hier, molecule=obj,
            resolution=self.resolution).get_selected_particles()
        return particles

    def create_membrane_density(self, file_out='membrane_localization.mrc'):
        """Create an MRC density file to visualize the membrane."""
        offset = 5.0 * self.thickness
        apix = 3.0
        resolution = 5.0

        # Create a density header of the requested size
        bbox = IMP.algebra.BoundingBox3D(
            IMP.algebra.Vector3D(-self.center - offset, -self.center - offset,
                                 -self.center - offset),
            IMP.algebra.Vector3D(self.center + offset, self.center + offset,
                                 self.center + offset))
        dheader = IMP.em.create_density_header(bbox, apix)
        dheader.set_resolution(resolution)
        dmap = IMP.em.SampledDensityMap(dheader)

        for vox in range(dmap.get_header().get_number_of_voxels()):
            c = dmap.get_location_by_voxel(vox)
            if self._is_membrane(c[2]) == 1:
                dmap.set_value(c[0], c[1], c[2], 1.0)
            else:
                dmap.set_value(c[0], c[1], c[2], 0.0)

        IMP.em.write_map(dmap, file_out)

    def _is_membrane(self, z):
        if ((z-self.center) < self.thickness/2.0 and
                (z-self.center) >= -self.thickness/2.0):
            return 1
        else:
            return 0


class ResidueProteinProximityRestraint(IMP.pmi.restraints.RestraintBase):
    """Restrain residue/residues to bind to unknown location in a target"""

    def __init__(self, hier, selection, cutoff=6., sigma=3., xi=0.01,
                 resolution=1.0, weight=1.0, label=None):
        """
        Constructor
        @param hier        Hierarchy of the system
        @param selection   Selection of residues and target;
                           syntax is (prot, r1, r2, target_prot) or
                           (prot1, r1, r2, target_prot, target_r1, target_r2)
        @param cutoff      Distance cutoff between selected segment and target
                           protein
        @param sigma       Distance variance between selected fragments
        @param xi          Slope of a distance-linear scoring function that
                           funnels the score when the particles are too
                           far away
        @param resolution  Resolution at which to apply restraint
        @param weight      Weight of the restraint
        @param label       Extra text to label the restraint so that it is
                           searchable in the output
        """
        self.hier = hier
        m = self.hier.get_model()

        super().__init__(
            m, name="ResidueProteinProximityRestraint", label=label,
            weight=weight)

        self.cutoff = cutoff
        self.sigma = sigma
        self.xi = xi
        self.resolution = resolution

        # Check selection
        print('selection', selection, isinstance(selection, tuple))
        if (not isinstance(selection, tuple)
                and not isinstance(selection, list)):
            raise ValueError("Selection should be a tuple or list")
        if len(selection) < 4:
            raise ValueError(
                "Selection should be (prot, r1, r2, target_prot) or "
                "(prot1, r1, r2, target_prot, target_r1, target_r2)")

        # Selection
        self.prot1 = selection[0]
        self.r1 = int(selection[1])
        self.r2 = int(selection[2])

        self.prot2 = selection[3]
        if len(selection) == 6:
            self.tr1 = int(selection[4])
            self.tr2 = int(selection[5])

        if self.r1 == self.r2:
            sel_resi = IMP.atom.Selection(
                self.hier, molecule=self.prot1, residue_index=self.r1,
                resolution=self.resolution).get_selected_particles()
        else:
            sel_resi = IMP.atom.Selection(
                self.hier, molecule=self.prot1,
                residue_indexes=range(self.r1, self.r2+1, 1),
                resolution=self.resolution).get_selected_particles()

        if len(selection) == 4:
            sel_target = IMP.atom.Selection(
                self.hier, molecule=self.prot2,
                resolution=self.resolution).get_selected_particles()

        elif len(selection) == 6:
            sel_target = IMP.atom.Selection(
                self.hier, molecule=self.prot2,
                residue_indexes=range(self.tr1, self.tr2+1, 1),
                resolution=self.resolution).get_selected_particles()

        self.included_ps = sel_resi + sel_target

        # Setup restraint
        distance = 0.0
        slack = cutoff*2

        br = IMP.isd.ResidueProteinProximityRestraint(
            m, self.cutoff, self.sigma, self.xi, True,
            'ResidueProteinProximityRestraint')

        print('Selected fragment and target lengths:', len(sel_resi),
              len(sel_target))

        # Setup close pair container
        # Find close pair within included_resi and included_target
        lsa_target = IMP.container.ListSingletonContainer(m)
        lsa_target.add(IMP.get_indexes(sel_target))

        lsa_resi = IMP.container.ListSingletonContainer(m)
        lsa_resi.add(IMP.get_indexes(sel_resi))

        self.cpc = IMP.container.CloseBipartitePairContainer(lsa_resi,
                                                             lsa_target,
                                                             distance,
                                                             slack)

        br.add_pairs_container(self.cpc)

        br.add_contribution_particles(sel_resi, sel_target)

        # Compute interpolation parameters
        yi = ((cutoff**2/(2*sigma**2)
               - math.log(1/math.sqrt(2*math.pi*sigma*sigma))+cutoff*xi/2.)
              / (cutoff/2.))
        interpolation_factor = -(cutoff/2.)*(xi-yi)
        max_p = (math.exp(-((distance+slack)**2)/(2*sigma**2))
                 / math.sqrt(2*math.pi*sigma*sigma))
        max_score = -math.log(max_p)

        # Add interpolation parameters
        br.set_yi(yi)
        br.set_interpolation_factor(interpolation_factor)
        br.set_max_score(max_score)

        self.rs.add_restraint(br)

        self.restraint_sets = [self.rs] + self.restraint_sets[1:]

    def get_container_pairs(self):
        """ Get particles in the close pair container """
        return self.cpc.get_indexes()

    def get_output(self):
        output = {}
        score = self.weight * self.rs.unprotected_evaluate(None)
        output["ResidueProteinProximityRestraint_score_" + self.label] \
            = str(score)

        return output

class PMFRestraint(IMP.Restraint):
    """Restrain particles based on a PMF lookup table with dynamic pair updates"""

    def __init__(self, m, particles, pmf_file, distance_cutoff=None, weight=1.0, label="None"):
        IMP.Restraint.__init__(self, m, "PMFRestraint %1%")

        self.m = m
        self.weight = weight
        self.label = label

        # Create the C++ restraint
        self.pmf = IMP.isd.PMFRestraint(self.m, pmf_file, distance_cutoff or 0.0)

        if distance_cutoff is not None:
            # Set up containers for dynamic pair updates
            self.lsc = IMP.container.ListSingletonContainer(
                self.m, [p.get_index() for p in particles])
            self.cpc = IMP.container.ClosePairContainer(self.lsc, distance_cutoff)
            self.pmf.set_containers(self.lsc, self.cpc)
            print(f"Set up dynamic pair generation with distance cutoff {distance_cutoff}Å")
        else:
            # Manual pair addition
            if isinstance(particles[0], (list, tuple)):
                for pair in particles:
                    if len(pair) == 2:
                        self.pmf.add_particle_pair(pair[0], pair[1])
            else:
                for i in range(len(particles)):
                    for j in range(i + 1, len(particles)):
                        print(f"Adding pair: {particles[i].get_name()} - {particles[j].get_name()}")
                        self.pmf.add_particle_pair(particles[i], particles[j])

    def unprotected_evaluate(self, da):
        return self.weight * self.pmf.unprotected_evaluate(da)

    def do_get_inputs(self):
        return self.pmf.get_inputs()

    def add_to_model(self):
        IMP.pmi.tools.add_restraint_to_model(self.m, self)

    def get_output(self):
        output = {}
        score = self.weight * self.unprotected_evaluate(None)
        output["_TotalScore"] = str(score)
        output["PMFRestraint_" + self.label] = str(score)
        return output

import IMP.em
import numpy as np
class EMIlanRestraint(IMP.pmi.restraints.RestraintBase):
    """EM restraint using IMP.em BayesEM3D functions"""

    def __init__(self, hier, em_map_file, resolution=20.0, normalize_target=True,
                 sigma=1.0, weight=1.0, label=None):
        """
        Setup EM restraint using BayesEM3D functions with validation.
        """
        m = hier.get_model()
        super().__init__(m, name="EMIlanRestraint", label=label, weight=weight)

        self.hier = hier
        self.resolution = resolution
        self.sigma = sigma

        # Load experimental map with validation
        if isinstance(em_map_file, str):
            self.exp_map = IMP.em.read_map(em_map_file)
            print(f"Loaded EM map from: {em_map_file}")
        else:
            self.exp_map = em_map_file
            print("Using provided DensityMap object")

        # CRITICAL VALIDATION: Check if map is valid
        self._validate_experimental_map()

        # CRITICAL: Calculate RMS for proper correlation calculation
        self.exp_map.calcRMS()

        # Collect particles with proper XYZR and Mass setup
        self.particles = []
        self._setup_particles()

        if len(self.particles) == 0:
            raise ValueError("No XYZR particles found in hierarchy!")

        print(f"EMIlanRestraint: {len(self.particles)} particles, resolution={resolution}Å")

        # Validate particle positions vs map bounds
        self._validate_particle_positions()

        # Optional: Normalize target map intensities using BayesEM3D
        if normalize_target:
            print("Normalizing target map intensities using BayesEM3D...")
            self._normalize_target_intensities()

        # Validate setup before creating restraint
        self._run_initial_validation()

        # Create the actual scoring restraint
        em_scorer = EMIlanScoringFunction(self)
        self.rs.add_restraint(em_scorer)

        # Store for output/debugging
        self.last_ccc = 0.0
        self.last_score = 0.0

    def _validate_experimental_map(self):
        """Validate that the experimental map is valid"""
        print("=== Validating Experimental Map ===")

        header = self.exp_map.get_header()
        data = self.exp_map.get_data()

        # Check dimensions
        nx, ny, nz = header.get_nx(), header.get_ny(), header.get_nz()
        nvox = header.get_number_of_voxels()
        print(f"Map dimensions: {nx} x {ny} x {nz} = {nvox} voxels")

        if nvox == 0:
            raise ValueError("Experimental map has zero voxels!")

        # Check data validity
        data_min = np.min(data)
        data_max = np.max(data)
        data_mean = np.mean(data)
        data_std = np.std(data)

        print(f"Map data range: [{data_min:.6f}, {data_max:.6f}]")
        print(f"Map data mean: {data_mean:.6f}, std: {data_std:.6f}")

        if np.isnan(data_mean) or np.isinf(data_mean):
            raise ValueError("Experimental map contains NaN or Inf values!")

        if data_max == data_min:
            raise ValueError("Experimental map has constant values (no variation)!")

        # Check voxel size
        voxel_size = self.exp_map.get_spacing()
        print(f"Voxel size: {voxel_size:.3f} Å")

        if voxel_size <= 0:
            raise ValueError("Invalid voxel size!")

    def _setup_particles(self):
        """Setup particles with proper Mass and XYZR"""
        sel = IMP.atom.Selection(self.hier, resolution=IMP.atom.ALL_RESOLUTIONS)
        total_particles = 0

        for p in sel.get_selected_particles():
            total_particles += 1
            if IMP.core.XYZR.get_is_setup(p):
                # Ensure Mass is set up (required for BayesEM3D)
                if not IMP.atom.Mass.get_is_setup(p):
                    # Set mass based on radius (rough approximation)
                    radius = IMP.core.XYZR(p).get_radius()
                    estimated_mass = (radius**3)  # Volume-based
                    IMP.atom.Mass.setup_particle(p, estimated_mass)

                # Validate mass
                mass = IMP.atom.Mass(p).get_mass()
                if mass <= 0 or np.isnan(mass) or np.isinf(mass):
                    print(f"Warning: Invalid mass {mass} for particle {p.get_name()}, setting to 1.0")
                    IMP.atom.Mass(p).set_mass(1.0)

                self.particles.append(p)

        print(f"Particle setup: {len(self.particles)}/{total_particles} particles have XYZR")

    def _validate_particle_positions(self):
        """Check if particles are within reasonable bounds of the map"""
        print("=== Validating Particle Positions ===")

        # Get map bounding box
        bbox = IMP.em.get_bounding_box(self.exp_map)
        map_min = bbox.get_corner(0)
        map_max = bbox.get_corner(1)

        print(f"Map bounding box: [{map_min[0]:.1f}, {map_min[1]:.1f}, {map_min[2]:.1f}] to [{map_max[0]:.1f}, {map_max[1]:.1f}, {map_max[2]:.1f}]")

        # Check particle positions
        coords = []
        radii = []
        masses = []

        for p in self.particles:
            xyz = IMP.core.XYZ(p)
            xyzr = IMP.core.XYZR(p)
            mass = IMP.atom.Mass(p)

            coord = [xyz.get_x(), xyz.get_y(), xyz.get_z()]
            coords.append(coord)
            radii.append(xyzr.get_radius())
            masses.append(mass.get_mass())

        coords = np.array(coords)
        radii = np.array(radii)
        masses = np.array(masses)

        # Particle statistics
        particle_min = np.min(coords, axis=0)
        particle_max = np.max(coords, axis=0)
        particle_center = np.mean(coords, axis=0)

        print(f"Particle bounding box: [{particle_min[0]:.1f}, {particle_min[1]:.1f}, {particle_min[2]:.1f}] to [{particle_max[0]:.1f}, {particle_max[1]:.1f}, {particle_max[2]:.1f}]")
        print(f"Particle center: [{particle_center[0]:.1f}, {particle_center[1]:.1f}, {particle_center[2]:.1f}]")
        print(f"Radius range: [{np.min(radii):.1f}, {np.max(radii):.1f}] Å")
        print(f"Mass range: [{np.min(masses):.1f}, {np.max(masses):.1f}]")

        # Check overlap
        map_center = (np.array(map_min) + np.array(map_max)) / 2
        print(f"Map center: [{map_center[0]:.1f}, {map_center[1]:.1f}, {map_center[2]:.1f}]")

        # Distance between centers
        center_distance = np.linalg.norm(particle_center - map_center)
        map_size = np.linalg.norm(np.array(map_max) - np.array(map_min))
        print(f"Distance between centers: {center_distance:.1f} Å (map size: {map_size:.1f} Å)")

        if center_distance > map_size:
            print("Warning: Particles are very far from map center!")

    def _normalize_target_intensities(self):
        """Normalize target map intensities using BayesEM3D"""
        # Generate test model map first to check compatibility
        test_model_map = IMP.em.bayesem3d_get_density_from_particle(
            self.exp_map, self.particles, self.resolution
        )

        test_data = test_model_map.get_data()
        test_sum = np.sum(test_data)
        test_max = np.max(test_data)

        print(f"Test model map before normalization: sum={test_sum:.6f}, max={test_max:.6f}")

        if test_sum > 1e-10:  # Only normalize if model map has reasonable density
            IMP.em.bayesem3d_get_normalized_intensities(
                self.exp_map, self.particles, self.resolution
            )
            print("Target map normalization complete")
        else:
            print("Warning: Model map density too low, skipping normalization")

    def _run_initial_validation(self):
        """Run initial validation tests similar to the IMP test cases"""
        print("=== Running Initial Validation ===")

        # Test 1: Target map self-correlation
        ccc_self = IMP.em.bayesem3d_get_cross_correlation_coefficient(
            self.exp_map, self.exp_map
        )
        print(f"Target map self-correlation: {ccc_self:.8f}")

        if np.isnan(ccc_self) or ccc_self < 0.99:
            print(f"Warning: Target map self-correlation is {ccc_self}, expected ~1.0")

        # Test 2: Generate initial model map
        print("Generating initial model map...")
        initial_model_map = IMP.em.bayesem3d_get_density_from_particle(
            self.exp_map, self.particles, self.resolution
        )

        # Check if model map has any density
        model_data = initial_model_map.get_data()
        model_sum = np.sum(model_data)
        model_max = np.max(model_data)

        print(f"Initial model map: sum={model_sum:.6f}, max={model_max:.6f}")

        if model_sum == 0 or np.isnan(model_sum):
            raise ValueError("Initial model map is empty or contains NaN!")

        # Test 3: Model map self-correlation
        ccc_model_self = IMP.em.bayesem3d_get_cross_correlation_coefficient(
            initial_model_map, initial_model_map
        )
        print(f"Model map self-correlation: {ccc_model_self:.8f}")

        # Test 4: Initial cross-correlation
        initial_ccc = IMP.em.bayesem3d_get_cross_correlation_coefficient(
            self.exp_map, initial_model_map
        )
        print(f"Initial cross-correlation: {initial_ccc:.6f}")

        if np.isnan(initial_ccc):
            raise ValueError("Initial cross-correlation is NaN!")

        print("Initial validation passed")


class EMIlanScoringFunction(IMP.Restraint):
    """Scoring function following IMP.pmi patterns"""

    def __init__(self, em_restraint):
        IMP.Restraint.__init__(self, em_restraint.model, "EMIlanScoringFunction")
        self.em_restraint = em_restraint
        self._evaluation_count = 0
        self.model_map = None

        # Add caching to reduce redundant evaluations
        self._last_positions_hash = None
        self._cached_score = None
        self._cached_ccc = None

    def _get_positions_hash(self):
        """Create a hash of current particle positions for caching"""
        coords = []
        for p in self.em_restraint.particles:
            xyz = IMP.core.XYZ(p)
            coords.extend([xyz.get_x(), xyz.get_y(), xyz.get_z()])
        return hash(tuple(np.round(coords, 6)))  # Round to avoid floating point issues

    def unprotected_evaluate(self, da):
        """Evaluate with caching to avoid redundant calculations"""

        self._evaluation_count += 1

        # CHECK CACHE FIRST
        current_hash = self._get_positions_hash()
        if (self._last_positions_hash is not None and
            current_hash == self._last_positions_hash and
            self._cached_score is not None):

            # Return cached result - suppress debug output for cached results
            return self._cached_score

        # Generate model density
        self.model_map = IMP.em.bayesem3d_get_density_from_particle(
            self.em_restraint.exp_map,
            self.em_restraint.particles,
            self.em_restraint.resolution,
            1.0
        )

        # Quick validation of model map
        model_data = self.model_map.get_data()
        model_sum = np.sum(model_data)

        if model_sum == 0 or np.isnan(model_sum) or np.isinf(model_sum):
            print(f"Warning: Invalid model map at evaluation {self._evaluation_count}, sum={model_sum}")
            score = 1000.0
            ccc = 0.0
        else:
            # Calculate cross-correlation
            ccc = IMP.em.bayesem3d_get_cross_correlation_coefficient(
                self.em_restraint.exp_map,
                self.model_map
            )

            # Validate CCC
            if np.isnan(ccc) or np.isinf(ccc):
                print(f"Warning: Invalid CCC at evaluation {self._evaluation_count}, ccc={ccc}")
                score = 1000.0
                ccc = 0.0
            else:
                # Convert to score (ensure it's valid)
                score = 1000.0 * (1.0 - ccc)

                if np.isnan(score) or np.isinf(score):
                    print(f"Warning: Invalid score at evaluation {self._evaluation_count}, score={score}")
                    score = 1000.0

        # CACHE THE RESULTS
        self._last_positions_hash = current_hash
        self._cached_score = score
        self._cached_ccc = ccc

        # Store results
        self.em_restraint.last_ccc = ccc
        self.em_restraint.last_score = score

        # Reduce debug output frequency
        if self._evaluation_count % 500 == 1:  # Only every 500 evaluations
            print(f"EM Evaluation #{self._evaluation_count}: CCC = {ccc:.6f}, Score = {score:.2f}")

        return score

    def do_get_inputs(self):
        return self.em_restraint.particles

    def get_output(self):
        """Return output dict following IMP.pmi pattern like other restraints"""
        output = {}
        score = self.em_restraint.weight * self.unprotected_evaluate(None)
        output["_TotalScore"] = str(score)
        output["EMIlanRestraint_" + (self.em_restraint.label or "NoLabel")] = str(score)
        output["EMIlanRestraint_CCC_" + (self.em_restraint.label or "NoLabel")] = str(self.em_restraint.last_ccc)
        return output

import math
import numpy as np
import scipy.ndimage
import IMP
import IMP.core
import IMP.atom
import IMP.em
import IMP.pmi.tools
import IMP.pmi.restraints

class EMArthurRestraint(IMP.Restraint):
    """EM restraint building a model density via weighted histogram + Gaussian blur
       and scoring with IMP.em.get_coarse_cc_coefficient.
    """

    def __init__(self, model, hierarchy, exp_map, resolution=20.0,
                 weight=1.0, label="EMArthur"):
        super().__init__(model, "EMArthurRestraint %1%")
        self.model = model
        self.hierarchy = hierarchy
        self.weight = weight
        self.label = label
        self.resolution = resolution
        self.threshold = 0.0

        if isinstance(exp_map, str):
            self.exp_map = IMP.em.read_map(exp_map)
        else:
            self.exp_map = exp_map
        self.exp_map.calcRMS()

        header = self.exp_map.get_header()
        self.voxel_size = header.get_spacing()
        self.nx = header.get_nx()
        self.ny = header.get_ny()
        self.nz = header.get_nz()

        sel = IMP.atom.Selection(hierarchy, resolution=IMP.atom.ALL_RESOLUTIONS)
        self.particles = []
        for p in sel.get_selected_particles():
            if IMP.core.XYZR.get_is_setup(p):
                if not IMP.atom.Mass.get_is_setup(p):
                    radius = IMP.core.XYZR(p).get_radius()
                    IMP.atom.Mass.setup_particle(p, radius**3)
                self.particles.append(p)
        if not self.particles:
            raise ValueError("EMArthurRestraint: no XYZR particles found in hierarchy")

        half_x = 0.5 * self.nx * self.voxel_size
        half_y = 0.5 * self.ny * self.voxel_size
        half_z = 0.5 * self.nz * self.voxel_size
        self.bins = (
            np.linspace(-half_x, half_x, self.nx + 1),
            np.linspace(-half_y, half_y, self.ny + 1),
            np.linspace(-half_z, half_z, self.nz + 1),
        )
        self.sigma = self.resolution / (4.0 * math.sqrt(2.0 * math.log(2.0))) / self.voxel_size

        self.coords = np.zeros((len(self.particles), 3), dtype=np.float32)
        self.weights = np.zeros(len(self.particles), dtype=np.float32)
        self.model_map = IMP.em.SampledDensityMap(header)
        self._eval_count = 0

    def _update_arrays(self):
        for i, p in enumerate(self.particles):
            xyz = IMP.core.XYZ(p)
            self.coords[i, 0] = xyz.get_x()
            self.coords[i, 1] = xyz.get_y()
            self.coords[i, 2] = xyz.get_z()
            self.weights[i] = IMP.atom.Mass(p).get_mass()

    def _fill_model_map(self, blurred):
        # get_data() returns 3D array matching map dimensions
        data_ptr = self.model_map.get_data()

        # Check if we can assign directly (shapes must match)
        if data_ptr.shape == blurred.shape:
            data_ptr[:] = blurred.astype(np.float64)
        else:
            # Reshape blurred to match data_ptr's shape
            # IMP might expect (nz, ny, nx) or (nx, ny, nz)
            data_ptr[:] = blurred.T.astype(np.float64)  # Try transpose first

        self.model_map.calcRMS()

    def unprotected_evaluate(self, da):
        self._eval_count += 1

        self._update_arrays()
        hist, _ = np.histogramdd(self.coords, bins=self.bins, weights=self.weights)
        hist = np.swapaxes(hist, 0, 2)
        blurred = scipy.ndimage.gaussian_filter(hist, self.sigma, truncate=4).astype(np.float32)

        if not np.isfinite(blurred.sum()):
            raise ValueError("EMArthurRestraint: invalid density values")
        self._fill_model_map(blurred)

        ccc = IMP.em.get_coarse_cc_coefficient(self.exp_map, self.model_map,
                                               self.threshold, False)
        if not np.isfinite(ccc):
            raise ValueError("EMArthurRestraint: invalid CCC value")
        score = self.weight * (1.0 - ccc)

        if self._eval_count % 500 == 1:
            print(f"EMArthur eval #{self._eval_count}: CCC={ccc:.6f}")
        return score

    def do_get_inputs(self):
        return self.particles

    def add_to_model(self):
        IMP.pmi.tools.add_restraint_to_model(self.model, self)

class EMArthurRestraintWrapper(IMP.pmi.restraints.RestraintBase):
    """PMI-friendly wrapper mirroring other basic restraints."""

    def __init__(self, model, hierarchy, exp_map, resolution=20.0,
                 weight=1.0, label=None):
        super().__init__(model, name="EMArthurRestraint",
                         label=label, weight=weight)
        restraint = EMArthurRestraint(model, hierarchy, exp_map,
                                      resolution=resolution,
                                      weight=weight,
                                      label=label or "EMArthur")
        self.rs.add_restraint(restraint)
        self._restraint = restraint

    def get_output(self):
        score = self._restraint.unprotected_evaluate(None)
        return {
            "_TotalScore": str(score),
            f"EMArthurRestraint_{self.label}": str(score)
        }
