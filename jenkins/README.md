ewoms jenkins build scripts:
--------------------------------

[B]build-ewoms.sh[/B]:
This is a helper script which contains functions for building,
testing and cloning ewoms and its dependencies.

[B]build.sh[/B]:
This script will build dependencies, then build ewoms and execute its tests.
It is intended for post-merge builds of the master branch.

[B]build-pr.sh[/B]:
This script will build dependencies, then build ewoms and execute its tests.
It inspects the $ghbPrBuildComment environmental variable to obtain a pull request
to use for opm-common, opm-parser, opm-material, opm-core and opm-grid (defaults to master)
and then builds $sha1 of ewoms.

It is intended for pre-merge builds of pull requests.

You can optionally specify a given pull request to use for opm-common, opm-parser, opm-material,
opm-core and dune-cornerpoint through the trigger.
The trigger line needs to contain opm-common=&lt;pull request number&gt; 
and/or opm-parser=&lt;pull request number&gt;
and/or opm-material=&lt;pull request number&gt; and/or opm-core=&lt;pull request number&gt;
and/or opm-grid=&lt;pull request number&gt;.
