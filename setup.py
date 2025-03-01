from setuptools import setup, find_packages

setup(
    name="ipa",
    version="0.1.0",
    author="Nikolai Bykov",
    author_email="nikolai.bykov@cnag.eu",
    description="Interaction Pattern Aggregation analysis (IPA)",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/encent/ipa",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "bioframe",
        "cooler",
        "cooltools",
        "pybbi",
        "tqdm",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    license="MIT",
    entry_points={
        'console_scripts': [
            'ipa=ipa.cli:main',
        ],
    },
)
