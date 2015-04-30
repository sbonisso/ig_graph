require 'rake/clean'  # for auto-cleaning
require 'open3'

namespace :iggraph do 
  num_proc = 8
  
  desc 'grid search over VJ/D k-mers for alleles'
  task :grid_search_k_allele do
    unless File.exists?("param_out_alleles.csv")
      data_f = "data/smab_data_mut10.fa"
      cmd = "./param_search.rb #{data_f} param_out_alleles.csv allele #{num_proc}"
      Open3.capture3(cmd)      
    end
  end

  desc 'plot VJ/D k-mer performance for alleles'
  task :plot_k_allele => :grid_search_k_allele do
    unless File.exists?("mat_allele.pdf")
      cmd = "./plot_matrix.R param_out_alleles.csv mat_allele.pdf heatmap"
      Open3.capture3(cmd)
    end
  end
  
  desc 'grid search over VJ/D k-mers for gene'
  task :grid_search_k_gene do
    unless File.exists?("param_out_gene.csv")
      data_f = "data/smab_data_mut10.fa"
      cmd = "./param_search.rb #{data_f} param_out_gene.csv gene #{num_proc}"
      Open3.capture3(cmd)
    end
  end

  desc 'plot VJ/D k-mer performance for genes'
  task :plot_k_gene => :grid_search_k_gene do
    unless File.exists?("mat_gene.pdf")
      cmd = "./plot_matrix.R param_out_gene.csv mat_gene.pdf heatmap"
      Open3.capture3(cmd)
    end
  end

end
