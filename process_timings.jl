using DataFrames

plot_data = false

tex_output = "pic/"
pic_output = "pictures/"

if plot_data
    using Gadfly
    set_default_plot_size(32cm, 18cm);

    width = 32cm
    height = 20cm

    pdf_width = 21cm
    pdf_height = 13cm
end

# Write the specified columns of a DataFrame to a csv file at the given
# location.
function save_data(data, columns::Array{Symbol,1}, filename::String)
    writetable("data/$filename.csv", data[columns]);
end

function save_data_by_algo(data::DataFrame, indices::Array{Symbol,1}, filename::String)
    save_data(data[data[:Algorithm] .== "mpz_probab_prime_p", :], indices, filename * "_gmp");
    save_data(data[data[:Algorithm] .== "Miller-Rabin", :], indices, filename * "_mr");
    save_data(data[data[:Algorithm] .== "Frobenius", :], indices, filename * "_frob");
end

function save_data_by_set(data::DataFrame, indices::Array{Symbol,1}, filename::String)
    save_data(data[data[:Set] .== "primes", :], indices, filename * "_primes");
    save_data(data[data[:Set] .== "composites", :], indices, filename * "_composites");
    save_data(data[data[:Set] .== "Mersenne numbers", :], indices, filename * "_mersenne_numbers");
    save_data(data[data[:Set] .== "Mersenne primes", :], indices, filename * "_mersenne_primes");
end

# Given a DataFrame containing the raw measurements, compute confidence
# intervals for the timings and somehow merge the other data.
function process_raw_data(df)
    μ = mean(df[:Time])
    σ = std(df[:Time])
    iters = sum(df[:Iterations])
    mults_per_iter = sum(df[:Multiplications]) / sum(df[:Iterations])
    number_of_pprimes = sum(df[:IsPrime])
    return DataFrame(MeanTime = μ, Min = μ - 4σ, Max = μ + 4σ,
                     Iterations = iters, Multiplications = mults_per_iter,
                     NumPrimes = number_of_pprimes)
end

# Generate a graph of the mean time needed by the various algorithms versus the
# length of the input.
function plot_timings(data::DataFrame, show_error_bars::Bool)
    default_labels = [Guide.XLabel("Eingabelänge in Bits"), Guide.YLabel("Zeit in Sekunden"), Guide.ColorKey("Algorithmen")];
    if show_error_bars
        return plot(data, x=:Bits, y=:MeanTime, color=:Algorithm, ymin=:Min, ymax=:Max,
                    Geom.errorbar, Geom.point, Scale.y_log10, Scale.x_log2, default_labels...)
    else
        return plot(data, x=:Bits, y=:MeanTime, color=:Algorithm,
                    Geom.point, Scale.y_log10, Scale.x_log2, default_labels...)
    end
end

# Save a given plot to disk.
function save_plot(plot, filename::String)
    #draw(PDF(string(pic_output, filename, ".pdf"), pdf_width, pdf_height), plot);
    draw(PNG(string(pic_output, filename, ".png"), width, height), plot);
end

function plot_with_and_without_errorbars(data::DataFrame, filename::String)
    save_data_by_algo(data, [:Bits, :MeanTime], filename);
    if plot_data
        save_plot(plot_timings(data, true), filename * "_error");
        save_plot(plot_timings(data, false), filename);
    end
end

# Compare the runtime of mpz_probab_prime_p and the RQFT with the Miller-Rabin
# test.
function plot_normalized(data::DataFrame, filename::String)
    save_data_by_algo(data, [:Bits, :NormalizedTime], filename);
    if plot_data
        p = plot(data, x=:Bits, y=:NormalizedTime, color=:Algorithm,
                 Geom.point, Scale.x_log2, Geom.hline, yintercept=[1,2],
                 Guide.XLabel("Eingabelänge in Bits"), Guide.YLabel("Relative Laufzeit (Miller-Rabin=1)"),
                 Guide.ColorKey("Algorithmus"));
        save_plot(p, filename);
    end
end

# Show, how long a single iteration of the algorithms takes without the
# pre-computation.
function plot_iteration_time(data::DataFrame, filename::String)
    save_data_by_set(data, [:Bits, :MeanTime], filename);
    if plot_data
        p = plot(data, x=:Bits, y=:MeanTime, color=:Algorithm, Geom.point,
                 Scale.x_log2, Scale.y_log10,
                 Guide.XLabel("Eingabelänge in Bits"), Guide.YLabel("Zeit in Sekunden"),
                 Guide.ColorKey("Algorithmus"));
        save_plot(p, filename);
    end
end

info("Read file with timing information");
data = readtable("timings_20140708.csv", separator=',');

info("Process data")
timings = by(data, [:Bits, :Algorithm, :Set, :Mode, :HammingWeight], process_raw_data);

function make_positive(x)
    x[x .< 0] = 1e-10;
    x
end

info("Split into preparation only and full test data");
timings_prep = timings[timings[:Mode] .== "prep", :];
timings_full = timings[timings[:Mode] .== "full", :];

# Calculate normalized time for full Algorithm
timings_doubled = join(timings_full, timings_full[timings_full[:Algorithm] .== "Miller-Rabin", :],
                       on=[:Set, :Bits, :Mode, :HammingWeight], kind=:outer);
delete!(timings_doubled, [:Algorithm_1, :Min_1, :Max_1, :Iterations_1, :Multiplications_1, :NumPrimes_1])
timings_normalized = by(timings_doubled, [:Bits, :Algorithm, :Set, :Mode],
                        df -> DataFrame(NormalizedTime = df[:MeanTime] ./ df[:MeanTime_1]));

# Calculate time per iteration
timings_doubled = join(timings_full, timings_prep, on=[:Set, :Bits, :Algorithm, :HammingWeight], kind=:inner)
delete!(timings_doubled, [:Min_1, :Max_1, :Iterations_1, :Multiplications_1, :NumPrimes_1])
timings_per_iteration = by(timings_doubled, [:Bits, :Algorithm, :Set, :Mode],
                           df -> DataFrame(MeanTime = make_positive(df[:MeanTime] .- df[:MeanTime_1])));

# Outer join to add a column (MeanTime_1) containing the mean runtime of the
# Miller-Rabin test.
timings_doubled = join(timings_per_iteration, timings_per_iteration[timings_per_iteration[:Algorithm] .== "Miller-Rabin", :],
                       on=[:Set, :Bits], kind=:outer);
timings_iter_normed = by(timings_doubled, [:Bits, :Algorithm, :Set],
                         df -> DataFrame(NormalizedTime = df[:MeanTime] ./ df[:MeanTime_1]));

info("Separate further according to the input set");
# Separate the timings according to the various input sets
prep_primes = timings_prep[timings_prep[:Set] .== "primes", :];
prep_composites = timings_prep[timings_prep[:Set] .== "composites", :];
prep_mersenne_numbers = timings_prep[timings_prep[:Set] .== "Mersenne numbers", :];
prep_mersenne_primes = timings_prep[timings_prep[:Set] .== "Mersenne primes", :];

full_primes = timings_full[timings_full[:Set] .== "primes", :];
full_composites = timings_full[timings_full[:Set] .== "composites", :];
full_mersenne_numbers = timings_full[timings_full[:Set] .== "Mersenne numbers", :];
full_mersenne_primes = timings_full[timings_full[:Set] .== "Mersenne primes", :];

norm_primes = timings_normalized[timings_normalized[:Set] .== "primes", :];
norm_composites = timings_normalized[timings_normalized[:Set] .== "composites", :];
norm_mersenne_numbers = timings_normalized[timings_normalized[:Set] .== "Mersenne numbers", :];
norm_mersenne_primes = timings_normalized[timings_normalized[:Set] .== "Mersenne primes", :];

iter_primes = timings_per_iteration[timings_per_iteration[:Set] .== "primes", :];
iter_composites = timings_per_iteration[timings_per_iteration[:Set] .== "composites", :];
iter_mersenne_numbers = timings_per_iteration[timings_per_iteration[:Set] .== "Mersenne numbers", :];
iter_mersenne_primes = timings_per_iteration[timings_per_iteration[:Set] .== "Mersenne primes", :];

iter_norm_primes = timings_iter_normed[timings_iter_normed[:Set] .== "primes", :];
iter_norm_composites = timings_iter_normed[timings_iter_normed[:Set] .== "composites", :];
iter_norm_mersenne_numbers = timings_iter_normed[timings_iter_normed[:Set] .== "Mersenne numbers", :];
iter_norm_mersenne_primes = timings_iter_normed[timings_iter_normed[:Set] .== "Mersenne primes", :];

info("Plot data");
# All kinds of timings in one plot
plot_with_and_without_errorbars(timings_full, "all");
plot_iteration_time(timings_per_iteration, "all_iter");
plot_normalized(timings_normalized, "all_normalized");
plot_normalized(timings_iter_normed, "all_iter_normalized");

# The runtime for just the pre-computation (mostly trial division) for each of the sets
plot_with_and_without_errorbars(prep_primes, "prep_primes")
plot_with_and_without_errorbars(prep_composites, "prep_composites")
plot_with_and_without_errorbars(prep_mersenne_numbers, "prep_mersenne_numbers")
plot_with_and_without_errorbars(prep_mersenne_primes, "prep_mersenne_primes")

# The runtime of the full algorithm
plot_with_and_without_errorbars(full_primes, "primes")
plot_with_and_without_errorbars(full_composites, "composites")
plot_with_and_without_errorbars(full_mersenne_numbers, "mersenne_numbers")
plot_with_and_without_errorbars(full_mersenne_primes, "mersenne_primes")

# The runtime per iteration
plot_iteration_time(iter_primes, "iter_primes");
plot_iteration_time(iter_composites, "iter_composites");
plot_iteration_time(iter_mersenne_numbers, "iter_mersenne_numbers");
plot_iteration_time(iter_mersenne_primes, "iter_mersenne_primes");

# The full runtime relative to Miller-Rabin
plot_normalized(norm_primes, "normalized_primes");
plot_normalized(norm_composites, "normalized_composites");
plot_normalized(norm_mersenne_numbers, "normalized_mersenne_numbers");
plot_normalized(norm_mersenne_primes, "normalized_mersenne_primes");

# The runtime per iteration
plot_normalized(iter_norm_primes, "iter_norm_primes");
plot_normalized(iter_norm_composites, "iter_norm_composites");
plot_normalized(iter_norm_mersenne_numbers, "iter_norm_mersenne_numbers");
plot_normalized(iter_norm_mersenne_primes, "iter_norm_mersenne_primes");

# How the input length corresponds to the number of multiplications needed
frob = timings_full[timings_full[:Multiplications] .!= 0, :];
save_data_by_set(frob, [:Bits, :Multiplications], "multiplications");
if plot_data
    p = plot(frob, x=:Bits, y=:Multiplications, color=:Set,
             Guide.XLabel("Eingabelänge in Bits"), Guide.YLabel("Anzahl der Multiplikationen"),
             Guide.ColorKey("Eingabemenge"),
             Scale.x_log2, Scale.y_log2);
    save_plot(p, "multiplications");
end

# How number of multiplications correlate with runtime
save_data_by_set(frob, [:Multiplications, :MeanTime, :Bits], "time_vs_multiplications");
if plot_data
    p = plot(frob, x=:Multiplications, y=:MeanTime, color=:Set,
             Guide.XLabel("Anzahl der Multiplikationen"), Guide.YLabel("Zeit in Sekunden"),
             Guide.ColorKey("Eingabemenge"),
             Scale.x_log2, Scale.y_log2);
    save_plot(p, "time_vs_multiplications");
end

info("Done");

# Numbers which were classified as "probably prime" by at least one iteration of a test.
possibly_prime = timings_full[timings_full[:NumPrimes] .!= 0, :];

# Hack that uses the fact that false positives are rare.  Thus, we have a false
# positive, if at least one test found the number tested to be composite.
# Should one number fool all iterations of all the tests, it would not show up.
# As the probability of that occurring are incredibly low, I will not bother
# doing this properly.
false_positives = possibly_prime[possibly_prime[:NumPrimes] .!= possibly_prime[:Iterations], :];

num_false_positives = size(false_positives[:NumPrimes], 1);

mersenne_composites = timings[timings[:Set] .== "Mersenne numbers", :];
mersenne_composites = mersenne_composites[mersenne_composites[:NumPrimes] .!= mersenne_composites[:Iterations], :];

# Count the number of iterations of the tests on numbers that are known to be
# composite.  This ignores the iterations 
num_iterations = sum(timings[timings[:Set] .== "composites", :][:Iterations]) + sum(mersenne_composites[:Iterations]);
file = open(string(tex_output, "false_positives.tex"), "w");
println(file, "Bei den insgesamt \$\\num{$num_iterations}\$ Iterationen der verschiedenen ",
        "Tests, die bei den zusammengesetzten Zahlen durchgeführt wurden, sind nur ",
        "\$\\num{$num_false_positives}\$ false positives aufgetreten.")
close(file);
