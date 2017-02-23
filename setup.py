from setuptools import setup

setup(# package information
      name="om10",
      version="0.7",
      author="Phil Marshall",
      author_email="dr.phil.marshall@gmail.com",
      description="A toolkit to use the OM10 catalog",
      long_description=open("README.md").read(),
      url="https://github.com/drphilmarshall/OM10",
      packages=['om10'],
      package_dir={'om10':'om10'},
      include_package_data=True,
      package_data={'om10': ['example_data/*']},
      classifiers=[
          "Development Status :: 4 - Beta",
          "License :: OSI Approved :: MIT License",
          "Intended Audience :: Developers",
          "Intended Audience :: Science/Research",
          "Operating System :: OS Independent",
          "Programming Language :: Python",
      ],
      install_requires=["numpy", "scipy", "matplotlib", "astropy", "sklearn"],
)
