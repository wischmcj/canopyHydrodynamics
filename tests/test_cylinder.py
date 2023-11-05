from __future__ import annotations

from canhydro.Cylinder import Cylinder


class TestCylinder:

    def setup_method(self, method):
        print(f"Setting up {method}")
        self.cylinder = Cylinder()

    def teardown_method(self, method):
        print(f"Tearing down up {method}")
        del self.cylinder

