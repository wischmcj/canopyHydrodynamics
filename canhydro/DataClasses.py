from __future__ import annotations

from collections import defaultdict
from collections.abc import Callable
from dataclasses import dataclass

from shapely.geometry import Point, Polygon
import networkx as nx


class MetaManager:
    "Returns a manager instance for an instance of the given object"

    def __get__(self, obj, objtype):
        if obj is None:
            return Manager(objtype)
        else:
            raise AttributeError(f"Manger isn't accessible via {objtype} instances")


class Manager:
    "Manages a client object"

    _store = defaultdict(list)

    def __init__(self, client):
        self._client = client
        self._client_name = f"{client.__module__}.{client.__qualname__}"

    def create(self, **kwargs):
        self._store[self._client_name].append(self._client(**kwargs))

    def all(self):
        return (obj for obj in self._store[self._client_name])

    def filter(self, a_lambda):
        if a_lambda.__code__.co_name != "<lambda>":
            raise ValueError("a lambda required")

        return (
            obj
            for obj in self._store[self._client_name]
            if eval(a_lambda.__code__, vars(obj).copy())
        )


class Model:
    "An example 'client'. Metamanager gives it an instance of a manager to manage it"

    cylinders = MetaManager()

    def __init__(self, **kwargs):
        if type(self) is Model:
            raise NotImplementedError

        class_attrs = self.__get_class_attributes(type(self))

        self.__init_instance(class_attrs, kwargs)

    def __get_class_attributes(self, cls):
        attrs = vars(cls)
        if "cylinders" in attrs:
            raise AttributeError(
                'class {} has an attribute named "cylinders" of type "{}"'.format(
                    type(self), type(attrs["cylinders"])
                )
            )
        attrs = {
            attr: obj
            for attr, obj in vars(cls).items()
            if not attr.startswith("_") and not isinstance(obj, Callable)
        }
        return attrs

    def __init_instance(self, attrs, kwargs_dict):
        for key, item in kwargs_dict.items():
            if key not in attrs:
                raise TypeError(f'Got an unexpected key word argument "{key}"')
            if isinstance(item, type(attrs[key])):
                setattr(self, key, item)
            else:
                raise TypeError(f"Expected type {type(attrs[key])}, got {type(item)}")


class CylinderList(Model):
    name = ""
    id = int()
    color = ""


if __name__ == "__main__":
    from pprint import pprint

    class CylinderList(Model):
        id = int()
        name = str()
        color = str()

        # def __init__(self,other):
        #     self.other = other

        def myFunc(self, **args):
            CylinderList.cylinders.create(**{"id": 1, "name": "brad", "color": "redd"})
            # CylinderList.cylinders.create(**{"id": 2, "name": "sylvia", "color": "blue"})
            # CylinderList.cylinders.create(**{"id": 3, "name": "paul", "color": "red"})
            # CylinderList.cylinders.create(**{"id": 4, "name": "brandon", "color": "yello"})
            # CylinderList.cylinders.create(**{"id": 5, "name": "martin", "color": "green"})

    # toby = CylinderList()
    # toby.myFunc()
    # pprint([vars(obj) for obj in Data.cylinders.filter(lambda: id == 1)])

    #causes the creation of cylinder list objects withing cylinderList.cylinders
    # class UserClass():
    #     myList = CylinderList()

    #     def addCyl(self,**args):
    #         breakpoint()
    #         self.myList.cylinders.create(**args)
    

# myUserClass = UserClass()
# breakpoint()
# myUserClass.addCyl(**{"id": 1, "name": "brrad", "color": "red"})
# breakpoint()
            
    
# self.arr = np.genfromtxt(file, delimiter=",", skip_header=True)[0:, :-1]
# cylinders = [self.create_cyl(row) for row in self.arr]

# # CylinderList.cylinders.create(**{"id": 1, "name": "brad", "color": "red"})
# # CylinderList.cylinders.create(**{"id": 2, "name": "sylvia", "color": "blue"})
# # CylinderList.cylinders.create(**{"id": 3, "name": "paul", "color": "red"})
# # CylinderList.cylinders.create(**{"id": 4, "name": "brandon", "color": "yello"})
# # CylinderList.cylinders.create(**{"id": 5, "name": "martin", "color": "green"})
# CylinderList.cylinders.create(**{"id": 6, "name": "annie", "color": "gray"})
# print([vars(obj) for obj in CylinderList.cylinders.filter(lambda: id == 1)])
# # print([vars(obj) for obj in CylinderList.cylinders.filter(lambda: 1 <= id <= 2)])
# # print([vars(obj) for obj in CylinderList.cylinders.filter(lambda: color == "blue")])
# breakpoint()


@dataclass
class Projection:
    plane: str
    polygon: Polygon()
    base_vector: list[int]
    anti_vector: list[int]
    angle: int()


    
@dataclass
class QsmGraph(nx.Graph):
    adj_dict:dict 

    def __post_init__(self, data=None, val=None, **attr):
        super(QsmGraph, self).__init__() #call graph constructor
        self.val = val

    def load_qsm():
        R = {}
        sid = self.df[' ID']
        pid = self.df[' parentID']
        sid.min()

        #QSM's are projected in different directions by swapping x,y and z values
        #however, for out edge calcualtions we need the typical orientation
        if self.projection == 'XZ':
            hypo= self.dy
            a_leg = self.dx
            b_leg = self.dz
        elif self.projection == 'YZ':
            hypo= self.dx
            a_leg = self.dy
            b_leg = self.dz
        else: # 'XY'
            hypo= self.dz
            a_leg = self.dy
            b_leg = self.dx

        attr = []
        gr = nx.Graph()
        dg = nx.DiGraph()
        for idx,curr_cyl_data in self.df.iterrows():
            #print a progress update once every 10 thousand or so cylinders
            # if np.random.uniform(0,1,1) < 0.0001:
            #     log.info(self.fileName + ':  wdgraph - adding edge  {} \n'.format(np.round((idx/len(self.df[0]))*100,decimals=1)))
            #     print('wdgraph - adding edges completed {} \n'.format(np.round((idx/len(self.df[0]))*100,decimals=1)))

            # Our first cylinder has ID 0 and parentID -1.
            # As a result cyilinder 0 is represented by and edge from node 0 to 1, and so on
            child_node = sid[idx] + 1
            par_node = pid[idx] + 1
            curr_cyl_id = idx
            par_cyl_id = np.where(sid == pid[idx])

            len_idx = curr_cyl_data[12]
            radius_idx = curr_cyl_data[9]
            poly_idx = self.pSVXY[curr_cyl_id]
            vector_idx = self.aV[curr_cyl_id]

            run = math.sqrt( a_leg[curr_cyl_id]**2 + b_leg[curr_cyl_id]**2)
            rise = hypo[curr_cyl_id]
            if run==0:
                slope_idx = 1 # straightDown e.g. is in flow
            else:
                slope_idx = rise/run

            sa_idx = 2*np.pi*radius_idx*(radius_idx +len_idx) - 2*np.pi*radius_idx*radius_idx
            pa_idx = poly_idx.area
            vol_idx = 2*np.pi*len_idx*radius_idx*radius_idx
            sa_to_vol_idx = sa_idx/vol_idx
            ang_idx = np.arctan(slope_idx)
            bo_idx = curr_cyl_data[20]

            gr.add_edge(child_node,par_node,length = len_idx,
                                radius =radius_idx,
                                aV = vector_idx,
                                poly = poly_idx,
                                inFlowGrade = slope_idx,
                                pa = pa_idx,
                                sa = sa_idx,
                                ang = ang_idx,
                                vol = vol_idx ,
                                sa_to_vol = sa_to_vol_idx,
                                bo = bo_idx
                            )
            if slope_idx<(0-(1/6)):
                dg.add_edge(child_node,par_node,length = len_idx,
                                radius = radius_idx,
                                aV = vector_idx,
                                poly = poly_idx,
                                inFlowGrade = slope_idx,
                                pa = pa_idx,
                                sa = sa_idx,
                                ang = ang_idx,
                                vol = vol_idx ,
                                sa_to_vol = sa_to_vol_idx,
                                bo = bo_idx
                            )
            else:
                dg.add_edge(par_node,child_node,length = len_idx,
                                radius = radius_idx,
                                aV = vector_idx,
                                poly = poly_idx,
                                inFlowGrade = slope_idx,
                                pa = pa_idx,
                                sa = sa_idx,
                                ang = ang_idx,
                                vol = ang_idx ,
                                sa_to_vol = sa_to_vol_idx,
                                bo = bo_idx
                            )
        self.graph = gr
        self.diGraph = dg
