# hera_op

`hera_op` provides a convenient and flexible framework for developing data
analysis pipelines for operating on HERA data. It facilitates "offline
processing", and is portable enough to operate on computer clusters with
batch submission systems or on local machines.

# How It Works

The `hera_op` package uses the `makeflow` system, which is a part of the
[Cooperative Computing Tools
package](https://github.com/cooperative-computing-lab/cctools) developed by the
[Cooperative Computing Lab](http://ccl.cse.nd.edu). The `hera_op` package
essentially converts a pipeline defined in a configuration file into a format
that can be parsed by `makeflow`. This process is also aware of aspects specific
to HERA data, such as the polarization features of the data, in order to build
an appropriate software pipeline. Once the `makeflow` instructions file has been
generated, the `makeflow` program itself is used to execute the steps in the
pipeline.

There are generally XX steps required to "build a pipeline":

1. Write *task scripts* that will be executed by `makeflow` for a given stage in
the pipeline. These scripts should generally be as atomic as possible, and
perform only a single logical component of a pipeline (though it may in turn
call several supporting scripts or commands).
2. Write a *configuration file* which defines the order of tasks to be
completed. This configuration file defines the logical flow of the pipeline, as
well as prerequisites for each task. It also allows for defining compute and
memory requirements, for systems that support resource management.
3. Use the provided *build_makeflow_from_config.py* script to build a `makeflow`
instruction file that specifies the pipeline tasks applied to the data files.
4. Use the provided *launch_makeflow_cluster.py* or *launch_makeflow_local.py*
to execute the pipeline in either a cluster/batch schedule environment, or on a
local machine, respectively.
5. (Optional) Use the provided *clean_up_makeflow.py* to clean up the work
directory for makeflow. This will remove the wrapper scripts and output files,
and generate a single log file for all jobs in the makeflow.

# Installation

As mentioned above, `hera_op` uses `makeflow` as the backing pipeline management
software. As such, `makeflow` must be installed. To install `makeflow` in your
home directory:
```
git clone https://github.com/cooperative-computing-lab/cctools.git
cd cctools
./configure --prefix=${HOME}/cctools
make clean
make install
export PATH=${PATH}:${HOME}/cctools/bin
```
For convenience, it is helpful to add the `export` statement to your `.bashrc`
file, so that the `makeflow` are always on your `PATH`.

To install the `hera_op` package, simply:
```
python setup.py install
```

# Task Scripts and Config Files

For documentation on building task scripts, see [the task scipts docs
page](add_link_here). For documentation on config files, see [the config file
docs page](add_link_here).

