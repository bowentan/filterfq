# filterfq

The program for cleaning/filtereing raw fastq files, providing the following functionalities.
    1. Automatically check quality system used in the raw fastq file(s)
    2. Filter out reads that have a number of 'N' bases, low average quality or a number of low quality base
    3. Convert quality system to specified system

## Getting Started

To use the program, you can fork a copy to your local machine by typing the following command in your work directory:

```
git clone https://github.com/bowentan/filterfq.git
```
or downloading the source codes or binary executables from the release page.

For those who download the source codes, use the following command to compile the program after unzipping it and entering the directory

```
./configure --prefix=/your/package/path
make
make install
```

## Running the tests

Explain how to run the automated tests for this system

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

1.0.0

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Bowen Tan** - *Initial work* - [bowentan](https://github.com/bowentan)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc
