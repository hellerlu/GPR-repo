from tools.input_cmd_funcs import command
import numpy as np

class materials():

    # predefined materials with ['name', rel permittivity, conductivity, rel permeability, magnetic loss]
    ballast =   ['ballast',  6,  0,    1, 0]
    dry_sand =  ['dry_sand', 7,  0,    1, 0]
    # wet_sand =  ['wet_sand', 15, 0,    1, 0]
    # wet_clay =  ['wet_clay', 25, 0.01, 1, 0]
    concrete =  ['concrete', 8,  0.01, 1, 0]
    dry_wood =  ['dry_wood', 2,  0.01, 1, 0]
    pss =       ['pss',      7,  0,    1, 0]
    apshalt =   ['asphalt',  8,  0.01, 1, 0]

    def get_matList(self):
        members = [attr for attr in dir(self) if not callable(getattr(self,attr)) and not attr.startswith("__") ]
        mat_list = []
        for mem in members:
            mat_list.append(getattr(self, mem))
        return mat_list

    def write_materials(self):
        mat_list = self.get_matList()
        f = ""
        for m in mat_list:
            f += command('material',m[1],m[2],m[3],m[4],m[0])
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
        if self.material == "wood":
            n_sl = int(np.floor((self.domain_size[0]-self.wood[1]-self.dist_dom_sleeper) / self.dist_sleepers )) + 1
        elif self.material == "steel":
            n_sl = int(np.floor((self.domain_size[0]-self.steel_1[1]-2*self.steel_2[1]-self.dist_dom_sleeper) / self.dist_sleepers )) + 1
        elif self.material == "concrete":
            n_sl = int(np.floor((self.domain_size[0]-self.concrete[1]-self.dist_dom_sleeper) / self.dist_sleepers )) + 1

        return n_sl

    def write(self, mat):
        f = ""
        for i in range(self.n_sleepers()):
            f += command('box',round(self.dist_dom_sleeper+i*self.dist_sleepers,3),round(self.top_height - mat[2],3),0,round(self.dist_dom_sleeper+i*self.dist_sleepers+mat[1],3),self.top_height,self.domain_size[2],mat[0])
        return f

    def write_sleepers(self):
        if self.material == "wood":
            return self.write(self.wood)
        elif self.material == "steel":
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
        f = ""
        f += command('box',0,self.start_height,round((self.domain_size[2]-self.rail_dist-self.steel_3[1])/2,3),self.domain_size[0],round(self.start_height+self.steel_3[2],3),round((self.domain_size[2]-self.rail_dist+self.steel_3[1])/2,3),self.steel_3[0])
        f += command('box',0,round(self.start_height+self.steel_3[2],3),round((self.domain_size[2]-self.rail_dist-self.steel_2[1])/2,3),self.domain_size[0],round(self.start_height+self.steel_3[2]+self.steel_2[2],3),round((self.domain_size[2]-self.rail_dist+self.steel_2[1])/2,3),self.steel_2[0])
        f += command('box',0,round(self.start_height+self.steel_3[2]+self.steel_2[2],3),round((self.domain_size[2]-self.rail_dist-self.steel_1[1])/2,3),self.domain_size[0],round(self.start_height+self.steel_3[2]+self.steel_2[2]+self.steel_1[2],3),round((self.domain_size[2]-self.rail_dist+self.steel_1[1])/2,3),self.steel_1[0])
        
        f += command('box',0,self.start_height,round((self.domain_size[2]+self.rail_dist-self.steel_3[1])/2,3),self.domain_size[0],round(self.start_height+self.steel_3[2],3),round((self.domain_size[2]+self.rail_dist+self.steel_3[1])/2,3),self.steel_3[0])
        f += command('box',0,round(self.start_height+self.steel_3[2],3),round((self.domain_size[2]+self.rail_dist-self.steel_2[1])/2,3),self.domain_size[0],round(self.start_height+self.steel_3[2]+self.steel_2[2],3),round((self.domain_size[2]+self.rail_dist+self.steel_2[1])/2,3),self.steel_2[0])
        f += command('box',0,round(self.start_height+self.steel_3[2]+self.steel_2[2],3),round((self.domain_size[2]+self.rail_dist-self.steel_1[1])/2,3),self.domain_size[0],round(self.start_height+self.steel_3[2]+self.steel_2[2]+self.steel_1[2],3),round((self.domain_size[2]+self.rail_dist+self.steel_1[1])/2,3),self.steel_1[0])
        
        return f