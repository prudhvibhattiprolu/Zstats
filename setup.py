import setuptools

with open("README.md", "r") as README:
    long_description = README.read()

setuptools.setup(
    name="Zstats",
    version="0.0.1",
    author="Prudhvi Bhattiprolu",
    author_email="prudhvibhattiprolu@gmail.com",
    description="Criteria for projected discovery and exclusion sensitivities of counting experiments",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    url="https://github.com/prudhvibhattiprolu",
    packages=setuptools.find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=["numpy>=1.19.1", "scipy>=1.4.1","mpmath>=1.1.0"],
)
