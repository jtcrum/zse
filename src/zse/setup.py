import setuptools

setuptools.setup(
    name="zeolite-simulation-environment-jcrum",  # Replace with your own username
    version="0.0.1",
    author="Jerry Crum",
    author_email="jcrum@nd.edu",
    description="A package with functions that automate zeolite strcuture building processes",
    url="https://github.com/jtcrum/zse",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
