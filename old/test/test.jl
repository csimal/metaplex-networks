using PkgTemplates

t = Template(;
        user="csimal",
        license="MIT",
        authors=["Cedric Simal"],
        plugins=[
            TravisCI(),
            Codecov(),
            Coveralls(),
            AppVeyor(),
            GitHubPages(),
           ],
       )

generate("TestPackage", t)
