using BinDeps
using Compat

@BinDeps.setup

libasl = library_dependency("libasl", aliases=["libasl.2", "libasl.2.0.3"])
libmp = library_dependency("libmp", aliases=["libmp.2", "libmp.2.0.3"])

# @osx_only begin
#   using Homebrew
#   provides(Homebrew.HB, "homebrew/science/asl", [libasl, libmp], os = :Darwin)
#   push!(Libdl.DL_LOAD_PATH, joinpath(Homebrew.prefix("asl"), "lib"))
# end

# Uncomment when there is a deb for the ASL.
# provides(AptGet, "libasl-dev", [libasl, libmp], os = :Linux)

@windows_only begin
  using WinRPM
  provides(WinRPM.RPM, "ampl-mp", [libasl, libmp], os = :Windows)
end

provides(Sources,
         URI("https://github.com/ampl/mp/archive/2.0.3.tar.gz"),
         [libasl, libmp],
         SHA="4ae38da883cfdf077d57c488b03756d9068b1d5b8552db983f6690246edc71a8",
         unpacked_dir="mp-2.0.3")

depsdir = BinDeps.depsdir(libasl)
prefix = joinpath(depsdir, "usr")
srcdir = joinpath(depsdir, "src", "mp-2.0.3")

provides(SimpleBuild,
         (@build_steps begin
            GetSources(libasl)
            (@build_steps begin
               ChangeDirectory(srcdir)
               (@build_steps begin
                  `wget https://gist.githubusercontent.com/dpo/dde4bf8030209fcf0569/raw/ed93e2653b51b5da754aabc89e08704421860009/a.diff`
                  `patch -p1 -i a.diff`
                  `wget https://github.com/ampl/mp/commit/ffede9ec6b131a3a8f8a35de9ba5bf4c648527b5.diff`
                  `patch -p1 -i ffede9ec6b131a3a8f8a35de9ba5bf4c648527b5.diff`
                  `cmake -DCMAKE_INSTALL_PREFIX=$prefix -DCMAKE_INSTALL_RPATH=$prefix/lib -DBUILD_SHARED_LIBS=True`
                  `make all`
                  `make test`
                  `make install`
                end)
             end)
          end), [libasl, libmp], os = :Unix)

@BinDeps.install @compat Dict(:libasl => :libasl)
