dncits <- c("RCOT", "CMIknn", "kpc_graph_2_10", "PCIT", "FCIT" )
settings <- data.frame(
    confounder = rep(c("age, sex", paste("age, sex, 10 PCs")), each=12),
    fct_relation = rep(rep(c("1", "2", "3"), each=4), 2),
    embedding_map = rep(c("same", "noisy", "noisy", "different"), 6),
    noise_x = rep(c(0, 50,500, 0), 6)
)

folder_path_T1E <- "C:\\Users\\Marco\\seadrive_root\\Marco Si\\Meine Bibliotheken\\Meine Bibliothek\\CITs\\P1 Deep CITs\\Paper_DNCITs\\Results\\CI\\rejection_rates"
sample_sizes <- c(30, 100, 300, 1000, 3000, 10000)
p_T1E <- create_test_plot(folder_path_T1E, dncits, settings, sample_sizes = sample_sizes)

folder_path_power <- "C:\\Users\\Marco\\seadrive_root\\Marco Si\\Meine Bibliotheken\\Meine Bibliothek\\CITs\\P1 Deep CITs\\Paper_DNCITs\\Results\\No_CI\\rejection_rates"
betas <- c(100.00, 75.00, 50.00, 25.00, 10.00, 7.50, 5.00, 2.50, 1.00, 0.75, 0.50, 0.25)
p_power <- create_test_plot(folder_path_power, dncits, settings, betas=betas)

