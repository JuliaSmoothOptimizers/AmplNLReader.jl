using BinDeps
using Compat

@BinDeps.setup

libasl = library_dependency("libasl", aliases=["libasl.2", "libasl.2.1.0"])
libmp = library_dependency("libmp", aliases=["libmp.2", "libmp.2.1.0"])

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
         URI("https://github.com/ampl/mp/archive/2.1.0.tar.gz"),
         [libasl, libmp],
         SHA="57d17db3e70e4a643c1b2141766a000b36057c2eeebd51964f30e2f8a56ee4d6",
         unpacked_dir="mp-2.1.0")

depsdir = BinDeps.depsdir(libasl)
prefix = joinpath(depsdir, "usr")
srcdir = joinpath(depsdir, "src", "mp-2.1.0")

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
