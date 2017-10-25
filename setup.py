from setuptools import find_packages, setup

setup(name="RetroSynthesis",
      version="1.0.0",
      description="A retrosynthetic tool that can identify enzyme/reaction pairs that when added \
       to a desired organism would produce a target chemical compound",
      author=["Leanne Whitmore","Corey M. Hudson"],
      author_email=["coreymhudson@gmail.com",
                    "lwhitmo@sandia.gov"
                    ],
      platforms=["linux",
                 "osx"
                 ],
      license="BSD 3 clause",
      url="https://github.com/sandialabs/BioRetroSynth",
      classifiers=['Development Status :: 3 - Alpha',
                   'Intended Audience :: Developers, AI Researchers, Bioinformaticists, bioengineers',
		   'Topic :: integer linear programming :: flux balance analysis',
   		   'License :: BSD 3 clause',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 2.6',
                   'Programming Language :: Python :: 2.7',
		   ],
      keywords='bioengineering, integer linear programming',
      test_suite="tests",
      packages=find_packages(),
      )
