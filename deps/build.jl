using BinDeps
using Compat

@BinDeps.setup

libasl = library_dependency("libasl", aliases=["libasl.3", "libasl.3.1.0"])
libmp = library_dependency("libmp", aliases=["libmp.3", "libmp.3.1.0"])

# Hopeless.
# @osx_only begin
#   using Homebrew
#   provides(Homebrew.HB, "homebrew/science/asl", [libasl, libmp], os = :Darwin)
#   push!(Libdl.DL_LOAD_PATH, joinpath(Homebrew.prefix("asl"), "lib"))
# end

# Uncomment when there is a deb for the ASL.
# provides(AptGet, "libasl-dev", [libasl, libmp], os = :Linux)

# Outdated!
@windows_only begin
  using WinRPM
  provides(WinRPM.RPM, "ampl-mp", [libasl, libmp], os = :Windows)
end

provides(Sources,
         URI("https://github.com/ampl/mp/archive/3.1.0.tar.gz"),
         [libasl, libmp],
         SHA="587c1a88f4c8f57bef95b58a8586956145417c8039f59b1758365ccc5a309ae9",
         unpacked_dir="mp-3.1.0")

depsdir = BinDeps.depsdir(libasl)
prefix = joinpath(depsdir, "usr")
srcdir = joinpath(depsdir, "src", "mp-3.1.0")

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

@BinDeps.install @compat Dict(:libasl => :libasl)
