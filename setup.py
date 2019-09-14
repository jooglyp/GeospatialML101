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
            "dask",
            "dask-ml",
            "scikit-learn",
            "dask_xgboost",
            "dask-searchcv",
            "fsspec",
            "imbalanced-learn",
            "seaborn",
            "graphviz",
            "pysal",
            "geopandas",
            "shapely",
            "geopy"
        ],
        entry_points={"console_scripts": ["geospatial_ml101=geospatial_ml101.main:main"]},
    )
