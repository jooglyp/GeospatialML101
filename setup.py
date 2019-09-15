import os

import setuptools

__folder__ = os.path.abspath(os.path.dirname(__file__))


if __name__ in ["__main__", "builtins"]:
    setuptools.setup(
        name="geospatial_ml101",
        include_package_data=True,
        packages=setuptools.find_packages(),
        python_requires=">=3.7",
        install_requires=[
            "numpy",
            "pandas",
            "fsspec",
            "seaborn",
            "graphviz",
            "pysal",
            "geopandas",
            "shapely",
            "geopy",
            "pyshp",
            "jupyterlab"
        ],
        entry_points={"console_scripts": ["data_processing=geospatial_ml101.main:data_processing",
                                          "data_generation=geospatial_ml101.main:data_generation"]},
    )
