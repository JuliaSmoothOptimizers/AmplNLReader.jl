using BinDeps

@BinDeps.setup

libasl = library_dependency("libasl", aliases=["libasl.2", "libasl.2.0.2"])
libmp = library_dependency("libmp", aliases=["libmp.2", "libmp.2.0.2"])

# Uncomment when the ASL makes it into Homebrew.jl.
# @osx_only begin
#   using Homebrew
#   provides(Homebrew.HB, "asl", [libasl, libmp], os = :Darwin)
# end

# Uncomment when there is a deb for the ASL.
# provides(AptGet, "libasl-dev", [libasl, libmp], os = :Linux)

@windows_only begin
  using WinRPM
  provides(WinRPM.RPM, "ampl-mp", [libasl, libmp], os = :Windows)
end

provides(Sources,
         URI("https://github.com/ampl/mp/archive/2.0.2.tar.gz"),
         [libasl, libmp],
         SHA="34d0ffef4af6f34ac5e50cfdb6819af823c2f6fb1db763b82de98b08d9069f63",
         unpacked_dir="mp-2.0.2")

depsdir = BinDeps.depsdir(libasl)
prefix = joinpath(depsdir, "usr")
srcdir = joinpath(depsdir, "src", "mp-2.0.2")

provides(SimpleBuild,
         (@build_steps begin
            GetSources(libasl)
            (@build_steps begin
               ChangeDirectory(srcdir)
               (@build_steps begin
                  `wget https://gist.githubusercontent.com/dpo/dde4bf8030209fcf0569/raw/ed93e2653b51b5da754aabc89e08704421860009/a.diff`
                  `patch -p1 -i a.diff`
                  `cmake -DCMAKE_INSTALL_PREFIX=$prefix -DCMAKE_INSTALL_RPATH=$prefix/lib -DBUILD_SHARED_LIBS=True`
                  `make all`
                  `make test`
                  `make install`
                end)
             end)
          end), [libasl, libmp], os = :Unix)

@BinDeps.install [:libasl => :libasl]
