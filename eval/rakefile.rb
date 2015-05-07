require 'rake/clean'  # for auto-cleaning
require 'open3'
require_relative 'run_data'

namespace :iggraph do 
  num_proc = 8
  
  smab_path = "#{ENV['HOME']}/git_repos/smab_lib/bin/smAb"
  data_f = "data/smab_data_mut10.fa"
  
  desc 'grid search over VJ/D k-mers for alleles'
  task :grid_search_k_allele do
    unless File.exists?("param_out_alleles.csv")
      cmd = "./param_search.rb #{data_f} param_out_alleles.csv allele #{num_proc}"
      Open3.capture3(cmd)      
    end
  end
  CLOBBER << "param_out_alleles.csv" # clobber calls clean

  desc 'plot VJ/D k-mer performance for alleles'
  task :plot_k_allele => :grid_search_k_allele do
    unless File.exists?("mat_allele.pdf")      
      check_cmd = "#{smab_path} check -f #{data_f}"
      out_r,err,pip = Open3.capture3(check_cmd)
      max_rate = out_r.chomp.to_f
      cmd = "./plot_matrix.R param_out_alleles.csv mat_allele.pdf heatmap #{max_rate}"
      Open3.capture3(cmd)
    end
  end
  CLEAN << "mat_allele.pdf"
  
  desc 'grid search over VJ/D k-mers for gene'
  task :grid_search_k_gene do
    unless File.exists?("param_out_gene.csv")
      cmd = "./param_search.rb #{data_f} param_out_gene.csv gene #{num_proc}"
      Open3.capture3(cmd)
    end
  end
  CLOBBER << "param_out_gene.csv"

  desc 'plot VJ/D k-mer performance for genes'
  task :plot_k_gene => :grid_search_k_gene do
    unless File.exists?("mat_gene.pdf")
      check_cmd = "#{smab_path} check -f #{data_f}"
      out_r,err,pip = Open3.capture3(check_cmd)
      max_rate = out_r.chomp.to_f
      cmd = "./plot_matrix.R param_out_gene.csv mat_gene.pdf heatmap #{max_rate}"
      Open3.capture3(cmd)
    end
  end
  CLEAN << "mat_gene.pdf"

  desc 'run on smAb data'
  task :test_smab_data do 
    unless File.exists?("smab_runs.csv")
      tools = ["igblast", "iggraph", "ihmmune"]
      File.open("smab_runs.csv", "w") do |f|
        run_data_supervised("./data/smab_data_mut10.fa", f, 
                            "./smab_run_pred", tools, "allele")
        #run_data_supervised("../tests/data/ten_ab.fa", f, "./smab_run_pred", "allele")
      end
    end
  end
  CLEAN << "smab_runs.csv"
  CLEAN << Dir.glob("smab_run_pred*")

  desc 'compare smAb data predictions'
  task :cmp_smab_data => :test_smab_data do 
    out_file_base = "smab_data_comparison"
    type = "allele"
    unless !Dir.glob(out_file_base+"*").empty?
      #
      pred_lst = Dir.glob("./smab_run_pred*")
      compare_unsupervised(pred_lst, out_file_base, "allele")
      
      # #
      # f_lst = Dir.glob("./smab_run_pred*")
      # tmp_f = Tempfile.new("filenames")
      # tmp_f.puts f_lst
      # tmp_f.close
      # #
      # name_f = Tempfile.new("names")
      # f_lst.each do |s| 
      #   name_f.puts s.split(".tab")[0].split("_")[-1].upcase        
      # end
      # name_f.close
      # FileUtils.cp(name_f.path, "/tmp/test_names.txt")
      # FileUtils.cp(tmp_f.path, "/tmp/test_file.txt")
      # #
      # valve_bin = "/home/stef/git_repos/valve/bin/valve"
      # valve_cmd = 
      #   "#{valve_bin} unsupervised -f #{tmp_f.path} -n #{name_f.path} -y #{type} -o #{out_file_base}"
      # Open3.capture3(valve_cmd)
    end
  end
  CLOBBER << Dir.glob("smab_data_comparison*")

  desc 'run on stanford_s22 data'
  task :test_stanford_s22_data do 
    unless File.exists?("stanford_s22_runs.csv")
      #tools = ["igblast", "iggraph", "ihmmune"]
      tools = ["igblast", "iggraph"]
      #tools = ["igblast"]
      File.open("stanford_s22_runs.csv", "w") do |f|
        run_data_unsupervised("./data/Stanford_S22_upcase.fasta",
                              f, 
                              "./stanford_s22_run_pred",
                              tools, 
                              "allele")
      end
    end
  end
  CLEAN << "stanford_s22_runs.csv"
  CLEAN << Dir.glob("stanford_s22_run_pred*")
  
  
  desc 'compare Stanford_S22 data predictions'
  task :cmp_stanford_s22_data => :test_stanford_s22_data do
    out_file_base = "stanford_s22_data_comparison"
    type = "allele"
    unless !Dir.glob(out_file_base+"*").empty?
      #
      puts 'hi'
      pred_lst = Dir.glob("./stanford_s22_run_pred*")
      compare_unsupervised(pred_lst, out_file_base, "allele")
    end
  end
  CLOBBER << Dir.glob("stanford_s22_data_comparison*")

  desc 'compare all Stanford_S22 tools'
  task :cmp_all_stanford_s22 do 
    out_file_base = "stanford_s22_all_comparison"
    name_re = Regexp.new(/(\w+)\_S22\_results\.txt$/)
    unless !Dir.glob(out_file_base + "*").empty?
      ["gene", "allele"].each do |y|
        out_file_base_s = [out_file_base, y].join("_")
        pred_lst = Dir.glob("./data/s22_results/*results.txt")        
        #pred_lst = Dir.glob("./data/s22_results_few/*results.txt")
        puts out_file_base_s
        
        name_lst = pred_lst.map{|s| name_re.match(s)[1]}
        pred_lst.push("./data/stanford_s22_run_pred_iggraph.tab")
        name_lst.push("IgGraph")
        puts pred_lst.to_s
        puts name_lst.to_s
        compare_unsupervised(pred_lst, 
                             out_file_base_s, 
                             y, 
                             name_lst)
      end
    end
  end

  desc 'plot all Stanford_S22'
  task :plot_all_stanford_s22 => :cmp_all_stanford_s22 do 
    unless File.exists?("stanford_s22_mat_V.pdf")
      cmd = "./plot_matrix.R stanford_s22_all_comparison_allele_V.csv stanford_s22_mat_V.pdf half_mat"
      Open3.capture3(cmd)
    end
  end
  

end
