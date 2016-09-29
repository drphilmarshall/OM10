from distutils.core import setup

setup(# package information
      name="om10",
      version="0.0.1dev",
      description='A toolkit to use the OM10 catalog ...',
      long_description=""" Not yet ... """,
      # What code to include as packages
      packages=['om10'],
      package_dir={'om10':'om10'},
      # What data to include as packages
      include_package_data=True,
      package_data={'om10': ['example_data/*']}
      )

