from tools.input_cmd_funcs import command
import numpy as np

class materials():

    def __init__(self, 
                 er_ballast = 6.5,
                 er_ballast_mix = 3.77, # USE dielectric_mixing.ipynb
                 er_dry_sand = 5,
                 pss_peplinski = [0.9, 0.1, 0.001, 0.002] # [sand fraction, clay fraction, lower water fraction, upper water fraction]
                 ):
        # Materials, add more in __init__ if more parameters should be varied
        self.ballast =       ['ballast',  er_ballast,  0,    1, 0]
        self.ballast_mix =   ['ballast_mix', er_ballast_mix , 0, 1, 0] 
        self.dry_sand =      ['dry_sand', er_dry_sand,  0,    1, 0]
        # Soil peplinski with ['name', sand fraction, clay fraction, bulk density g/cm3, density g/cm3 of sand, lower volumetry water fraction, upper vol water fraction]
        self.pss =           ['pss',  pss_peplinski[0], pss_peplinski[1], 2, 2.66, pss_peplinski[2], pss_peplinski[3]]


        # predefined materials with ['name', rel permittivity, conductivity, rel permeability, magnetic loss]
        # wet_sand =    ['wet_sand', 15, 0,    1, 0]
        # wet_clay =    ['wet_clay', 25, 0.01, 1, 0]
        self.concrete =      ['concrete', 8,  0.01, 1, 0]
        self.dry_wood =      ['dry_wood', 2,  0.01, 1, 0]
        self.asphalt =       ['asphalt',  8,  0.01, 1, 0]
        self.gravel =        ['gravel',   5,  0,    1, 0]
        self.fouling =       ['fouling',  5,  0,    1, 0]    
        self.fouling_mix =   ['fouling_mix',6.1, 0.01, 1, 0]

    
    ## Getting the material list automatically:
    # def get_matList(self):
    #     members = [attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__") ]
    #     mat_list = []
    #     for mem in members:
    #         mat_list.append(getattr(self, mem))
    #     return mat_list

    def write_materials(self):
        # If function get_matList is active:
        # mat_list = self.get_matList()

        mat_list = [[self.ballast,self.ballast_mix,self.dry_sand,self.concrete,self.dry_wood,self.asphalt,self.gravel,self.fouling, self.fouling_mix],[self.pss]]

        f = ""
        for m in mat_list[0]:
            f += command('material',m[1],m[2],m[3],m[4],m[0])
        for m in mat_list[1]:
            f += command('soil_peplinski',m[1],m[2],m[3],m[4],m[5],m[6],m[0])
        return f
            
class sleepers():

    def __init__(self, material, dist_dom_sleeper, dist_sleepers, top_height, domain_size):
            self.material = material
            self.dist_dom_sleeper = dist_dom_sleeper
            self.dist_sleepers = dist_sleepers
            self.top_height = top_height
            self.domain_size = domain_size

    # predefined sleepers with ['material', width [m], height[m]]
    wood =      ['dry_wood',0.15,   0.26]   # SBB Buchenschwelle RP IV K
    steel_1 =   ['pec',     0.220,  0.01]   # Österreichische Staatsbahnen, Platte oben approx
    steel_2 =   ['pec',     0.008,  0.1]    # Österreichische Staatsbahnen, Platten schräg seitlich approx
    concrete =  ['concrete',0.17,   0.1]

    def n_sleepers(self):
        """ Calculates number of sleepers in domain
    
        Input:
            self
        Output:
            n_sl:   Number of sleepers in domain
        """
        # Cases
        if self.material == "wood":
            n_sl = int(np.floor((self.domain_size[0]-self.wood[1]-self.dist_dom_sleeper) / self.dist_sleepers )) + 1
        elif self.material == "steel":
            n_sl = int(np.floor((self.domain_size[0]-self.steel_1[1]-2*self.steel_2[1]-self.dist_dom_sleeper) / self.dist_sleepers )) + 1
        elif self.material == "concrete":
            n_sl = int(np.floor((self.domain_size[0]-self.concrete[1]-self.dist_dom_sleeper) / self.dist_sleepers )) + 1

        return n_sl

    def write(self, mat):
        """ Writes sleeper command in an string array
    
        Input:
            self
        Output:
            f:   string containing all sleeper box commands
        """
        # Initialise    
        f = ""
        # Loop through number of sleepers
        for i in range(self.n_sleepers()):
            f += command('box',round(self.dist_dom_sleeper+i*self.dist_sleepers,3),round(self.top_height - mat[2],3),0,round(self.dist_dom_sleeper+i*self.dist_sleepers+mat[1],3),self.top_height,self.domain_size[2],mat[0])
        return f

    def write_sleepers(self):
        """ Writes sleeper commands depending on material to .in-file
    
        Input:
            self
        Output:
            f:    String containing all sleeper box commands depending on material
        """
        # Cases
        if self.material == "wood":
            return self.write(self.wood)
        elif self.material == "steel":
            # Steel is composed of three boxes
            f = ""
            f += self.write(self.steel_1)
            for i in range(self.n_sleepers()):
                f += command('box',round(self.dist_dom_sleeper+i*self.dist_sleepers- self.steel_2[1],3),round(self.top_height - self.steel_2[2],3),0,round(self.dist_dom_sleeper+i*self.dist_sleepers,3),self.top_height,self.domain_size[2],self.steel_2[0])
                f += command('box',round(self.dist_dom_sleeper+i*self.dist_sleepers+ self.steel_1[1],3),round(self.top_height - self.steel_2[2],3),0,round(self.dist_dom_sleeper+i*self.dist_sleepers+ self.steel_1[1] + self.steel_2[1],3),self.top_height,self.domain_size[2],self.steel_2[0])
            return f 
        elif self.material == "concrete":
            return self.write(self.concrete)

class rails():

    def __init__(self, top_height_sleepers, domain_size):
        self.start_height = top_height_sleepers
        self.domain_size = domain_size

    # Dimensions of rail according to SBB IV or 54 E2
    steel_1 =   ['pec', 0.067, 0.040]                   # top flange
    steel_2 =   ['pec', 0.016, 0.161 - 0.040 - 0.020]   # web
    steel_3 =   ['pec', 0.125, 0.020]                   # bottom flange

    rail_dist = 1.435 # Normalspur

    def write_rails(self):
        """ Writes rail command into string array
    
        Input:
            self
        Output:
            f:    String containing all rail box commands
        """

        f = ""
        f += command('box',0,self.start_height,round((self.domain_size[2]-self.rail_dist-self.steel_3[1])/2,3),self.domain_size[0],round(self.start_height+self.steel_3[2],3),round((self.domain_size[2]-self.rail_dist+self.steel_3[1])/2,3),self.steel_3[0])
        f += command('box',0,round(self.start_height+self.steel_3[2],3),round((self.domain_size[2]-self.rail_dist-self.steel_2[1])/2,3),self.domain_size[0],round(self.start_height+self.steel_3[2]+self.steel_2[2],3),round((self.domain_size[2]-self.rail_dist+self.steel_2[1])/2,3),self.steel_2[0])
        f += command('box',0,round(self.start_height+self.steel_3[2]+self.steel_2[2],3),round((self.domain_size[2]-self.rail_dist-self.steel_1[1])/2,3),self.domain_size[0],round(self.start_height+self.steel_3[2]+self.steel_2[2]+self.steel_1[2],3),round((self.domain_size[2]-self.rail_dist+self.steel_1[1])/2,3),self.steel_1[0])
        
        f += command('box',0,self.start_height,round((self.domain_size[2]+self.rail_dist-self.steel_3[1])/2,3),self.domain_size[0],round(self.start_height+self.steel_3[2],3),round((self.domain_size[2]+self.rail_dist+self.steel_3[1])/2,3),self.steel_3[0])
        f += command('box',0,round(self.start_height+self.steel_3[2],3),round((self.domain_size[2]+self.rail_dist-self.steel_2[1])/2,3),self.domain_size[0],round(self.start_height+self.steel_3[2]+self.steel_2[2],3),round((self.domain_size[2]+self.rail_dist+self.steel_2[1])/2,3),self.steel_2[0])
        f += command('box',0,round(self.start_height+self.steel_3[2]+self.steel_2[2],3),round((self.domain_size[2]+self.rail_dist-self.steel_1[1])/2,3),self.domain_size[0],round(self.start_height+self.steel_3[2]+self.steel_2[2]+self.steel_1[2],3),round((self.domain_size[2]+self.rail_dist+self.steel_1[1])/2,3),self.steel_1[0])
        
        return f