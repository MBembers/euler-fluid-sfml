# What does the program do?

This program is a visualisation of airflow in a wind tunnel.
The air is going from left to right. It carries a black smoke to show
the flow around a square box.

# How to use it?

You should have git and Cmake installed on your machine

1. If you use Linux, install SFML's dependencies using your system package manager. On Ubuntu and other Debian-based distributions you can use the following commands:
    ```
    sudo apt update
    sudo apt install \
        libxrandr-dev \
        libxcursor-dev \
        libudev-dev \
        libfreetype-dev \
        libopenal-dev \
        libflac-dev \
        libvorbis-dev \
        libgl1-mesa-dev \
        libegl1-mesa-dev
    ```
2. Build the project. Most popular IDEs support CMake projects with very little effort on your part.

    - [VS Code](https://code.visualstudio.com) via the [CMake extension](https://code.visualstudio.com/docs/cpp/cmake-linux)
    - [Visual Studio](https://docs.microsoft.com/en-us/cpp/build/cmake-projects-in-visual-studio?view=msvc-170)
    - [CLion](https://www.jetbrains.com/clion/features/cmake-support.html)

    Using CMake from the command line is straightforward as well.

    For a single-configuration generator (typically the case on Linux and macOS):

    ```
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build
    ```

    For a multi-configuration generator (typically the case on Windows):

    ```
    cmake -S . -B build
    cmake --build build --config Release
    ```

3. Executable will appear in ./build/bin directory
