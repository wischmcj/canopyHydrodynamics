from __future__ import annotations

from canhydro.global_vars import DIR, test_input_dir

DIR = DIR
test_input_dir = test_input_dir

# def test_create_cylinders(basic_forest):
#     actual = basic_forest.get_collection_data("1_TenCyls.csv")
#     expected = ten_cyls_rows
#     assert expected == actual


def test_sandpit(ten_cyls_col):
    ten_cyls_col.project_cylinders(plane="XZ")
    ten_cyls_col.initialize_graph()
    assert 1 == 1
