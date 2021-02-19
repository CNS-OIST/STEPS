node {
    library identifier: 'bbp@master', retriever: modernSCM(
        [$class:'GitSCMSource',
         remote: 'ssh://bbpcode.epfl.ch/hpc/jenkins-pipeline'])

    spack("steps",
          "git@github.com:CNS-OIST/HBP_STEPS.git",
          "+distmesh",
          test: "root")
}
