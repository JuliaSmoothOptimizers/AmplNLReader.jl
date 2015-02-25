using BinDeps

@BinDeps.setup

libasl = library_dependency("libasl", aliases=["libasl.1", "libasl.1.3.0"])
libmp = library_dependency("libmp", aliases=["libmp.1", "libmp.1.3.0"])

# Uncomment when the ASL makes it into Homebrew.jl.
# @osx_only begin
#     using Homebrew
#     provides(Homebrew.HB, "asl", [libasl, libmp], os = :Darwin)
# end

# Uncomment when there is a deb for the ASL.
# provides(AptGet, "libasl-dev", [libasl, libmp], os = :Linux)

# Uncomment when there is a Windows RPM for the ASL.
# @windows_only begin
#   using WinRPM
#   provides(WinRPM.RPM, "asl", [libasl, libmp], os = :Windows)
# end

provides(Sources,
         URI("https://github.com/ampl/mp/archive/1.3.0.tar.gz"),
         [libasl, libmp],
         SHA="aacac2c8f697e40d355457555774b6c1f71f36575a32a40fabef80f466d991b6",
         unpacked_dir="mp-1.3.0")

depsdir = BinDeps.depsdir(libasl)
prefix = joinpath(depsdir, "usr")
srcdir = joinpath(depsdir, "src", "mp-1.3.0")

provides(SimpleBuild,
         (@build_steps begin
            GetSources(libasl)
            (@build_steps begin
               ChangeDirectory(srcdir)
               (@build_steps begin
                  # Until the next release, we need to patch os-test.cc
                  `wget https://github.com/ampl/mp/commit/a9049249d579b6270826287c7162ed47443a1716.diff`
                  `cat a9049249d579b6270826287c7162ed47443a1716.diff` |> `patch -p1`
                  `cmake -DCMAKE_INSTALL_PREFIX=$prefix -DBUILD_SHARED_LIBS=True`
                  `make all`
                  `make test`
                  `make install`
                end)
             end)
          end), [libasl, libmp], os = :Unix)

@BinDeps.install [:libasl => :libasl]
