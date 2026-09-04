from setuptools import setup, find_packages

setup(
    name="bbetd",
    version="1.0.0",
    description="Block-Based Extended Tandem Duplication gene identification pipeline",
    author="zhixue",
    author_email="xuehzh95@foxmail.com",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "biopython>=1.78",
    ],
    entry_points={
        "console_scripts": [
            "bbetd = bbetd.cli:main",
        ],
    },
)