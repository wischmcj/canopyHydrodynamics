from __future__ import annotations

from collections import defaultdict
from collections.abc import Callable
from dataclasses import dataclass

from shapely.geometry import Point, Polygon


class MetaManager:
    "Returns a manager instance for an instance of the given object"

    def __get__(self, obj, objtype):
        if obj is None:
            return Manager(objtype)
        else:
            raise AttributeError(
                f"Manger isn't accessible via {objtype} instances"
            )


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
                raise TypeError(
                    f"Expected type {type(attrs[key])}, got {type(item)}"
                )


class CylinderList(Model):
    name = ''
    id = int()
    color = ''


if __name__ == "__main__":
    from pprint import pprint

    class Data(Model):
        asdf = ''
        adf = int()
        asdfs = ''

        # def __init__(self,other):
        #     self.other = other

        def myFunc(self):
            breakpoint()
            Data.cylinders.create(**{"id": 1, "name": "brad", "color": "red"})
            Data.cylinders.create(**{"id": 2, "name": "sylvia", "color": "blue"})
            Data.cylinders.create(**{"id": 3, "name": "paul", "color": "red"})
            Data.cylinders.create(**{"id": 4, "name": "brandon", "color": "yello"})
            Data.cylinders.create(**{"id": 5, "name": "martin", "color": "green"})

    toby = Data()
    # toby.myFunc()
    breakpoint()
    pprint([vars(obj) for obj in Data.cylinders.filter(lambda: id == 1)])
    breakpoint()

# CylinderList.cylinders.create(**{"id": 1, "name": "brad", "color": "red"})
# CylinderList.cylinders.create(**{"id": 2, "name": "sylvia", "color": "blue"})
# CylinderList.cylinders.create(**{"id": 3, "name": "paul", "color": "red"})
# CylinderList.cylinders.create(**{"id": 4, "name": "brandon", "color": "yello"})
# CylinderList.cylinders.create(**{"id": 5, "name": "martin", "color": "green"})
# CylinderList.cylinders.create(**{"id": 6, "name": "annie", "color": "gray"})
# print([vars(obj) for obj in CylinderList.cylinders.filter(lambda: id == 1)])
# print([vars(obj) for obj in CylinderList.cylinders.filter(lambda: 1 <= id <= 2)])
# print([vars(obj) for obj in CylinderList.cylinders.filter(lambda: color == "blue")])


@dataclass
class Projection:
    plane: str
    polygon: Polygon()
    base_vector: list[int]
    anti_vector: list[int]
    angle: int()