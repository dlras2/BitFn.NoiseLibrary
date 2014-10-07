using System.Diagnostics.Contracts;

namespace BitFn.NoiseLibrary
{
	public interface INoiseService4D
	{
		[Pure]
		double ComputeNoise(double x, double y, double z, double w);
	}
}