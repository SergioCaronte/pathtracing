#include "multijittered.h"

#define RAND_FLOAT(a) (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * a


MultiJittered::MultiJittered(void)
	: Sampler()
{}

MultiJittered::MultiJittered(const int num_samples)
	: 	Sampler(num_samples)
{
	generate_samples();
}

MultiJittered::MultiJittered(const int num_samples, const int m)
	: 	Sampler(num_samples, m)
{
	generate_samples();
}

MultiJittered::MultiJittered(const MultiJittered& mjs)
	: Sampler(mjs)
{
	generate_samples();
}


MultiJittered&
MultiJittered::operator= (const MultiJittered& rhs)
{
	if (this == &rhs)
		return (*this);

	Sampler::operator=(rhs);

	return (*this);
}

MultiJittered*
MultiJittered::clone(void) const
{
	return (new MultiJittered(*this));
}

MultiJittered::~MultiJittered(void) {}

void
MultiJittered::generate_samples(void)
{
	// num_samples needs to be a perfect square
	int n = (int)sqrt((float)num_samples);
	float subcell_width = 1.0 / ((float) num_samples);

	// fill the samples array with dummy points to allow us to use the [ ] notation when we set the
	// initial patterns

	Point fill_point;
	for (int j = 0; j < num_samples * num_sets; j++)
		samples.push_back(fill_point);

	// distribute points in the initial patterns

	for (int p = 0; p < num_sets; p++)
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++) {
				samples[i * n + j + p * num_samples].x = ((i * n + j) * subcell_width + RAND_FLOAT(subcell_width));
				samples[i * n + j + p * num_samples].y = ((j * n + i) * subcell_width + RAND_FLOAT(subcell_width));
			}

	// shuffle x coordinates

	for (int p = 0; p < num_sets; p++)
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				int k =  std::max(0, rand() % (n - j) + j - 1);
				float t = samples[i * n + j + p * num_samples].x;
				samples[i * n + j + p * num_samples].x = samples[i * n + k + p * num_samples].x;
				samples[i * n + k + p * num_samples].x = t;
			}

	// shuffle y coordinates

	for (int p = 0; p < num_sets; p++)
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				int k = std::max(0, rand() % (n - j) + j - 1);
				float t = samples[j * n + i + p * num_samples].y;
				samples[j * n + i + p * num_samples].y = samples[k * n + i + p * num_samples].y;
				samples[k * n + i + p * num_samples].y = t;
		}
}


