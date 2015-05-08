using BinDeps

@BinDeps.setup

libasl = library_dependency("libasl", aliases=["libasl.1", "libasl.1.4.0"])
libmp = library_dependency("libmp", aliases=["libmp.1", "libmp.1.4.0"])

# Uncomment when the ASL makes it into Homebrew.jl.
# @osx_only begin
#   using Homebrew
#   provides(Homebrew.HB, "asl", [libasl, libmp], os = :Darwin)
# end

# Uncomment when there is a deb for the ASL.
# provides(AptGet, "libasl-dev", [libasl, libmp], os = :Linux)

# Uncomment when there is a Windows RPM for the ASL.
# @windows_only begin
#   using WinRPM
#   provides(WinRPM.RPM, "asl", [libasl, libmp], os = :Windows)
# end

provides(Sources,
         URI("https://github.com/ampl/mp/archive/2.0.0.tar.gz"),
         [libasl, libmp],
         SHA="b8bc51dfbf3db1628e0eb029a3f1305b0640018f67ff6cc397bd06a1e63a1e09",
         unpacked_dir="mp-2.0.0")

depsdir = BinDeps.depsdir(libasl)
prefix = joinpath(depsdir, "usr")
srcdir = joinpath(depsdir, "src", "mp-2.0.0")

provides(SimpleBuild,
         (@build_steps begin
            GetSources(libasl)
            (@build_steps begin
               ChangeDirectory(srcdir)
               (@build_steps begin
                  `cmake -DCMAKE_INSTALL_PREFIX=$prefix -DCMAKE_INSTALL_RPATH=$prefix/lib -DBUILD_SHARED_LIBS=True`
                  `make all`
                  `make test`
                  `make install`
                end)
             end)
          end), [libasl, libmp], os = :Unix)

@BinDeps.install [:libasl => :libasl]
