/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2014-2017 Zuse Institute Berlin                             */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CHECK_POINT_HH
#define CHECK_POINT_HH

#include <string>
#include <sstream>
#include <memory>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <vector>
#include <utility>

#include <boost/filesystem.hpp>

#include <dune/grid/common/backuprestore.hh>
#include <dune/istl/bvector.hh>

#include <fem/gridmanager.hh>
#include <fem/variables.hh>

namespace Kaskade {
  struct CheckPoint
  {
    /**
     * @brief CheckPoint creates a new empty directory with current time stamp as name.
     *  Later writes will refer to this dirctory.
     */
    CheckPoint();

    /**
     * @brief CheckPoint. If the directory with pathName already exists later reads will refer to this directory.
     * Otherwise a new directory will be created where current time stamp is appended to pathName. This can then be used for writing.
     *
     * @param pathName is the directory which is to be created (if not existing) or where is to be read from (if existing).
     */
    CheckPoint(std::string const& pathName);

    /**
     * @brief readMode returns true if this check point was created in read mode (i.e. with existing directory given).
     */
    bool readMode() const {
      return timeStamp_.empty();
    }

    /**
     * @brief timeStamp gives the time stamp when directory was created (at construction of this CheckPoint instance), empty when directory already existed.
     * @return time stamp.
     */
    std::string timeStamp() const {
      return timeStamp_;
    }

    /**
     * @brief writes grid into file with name gridName.
     *
     * Does not work for UGGrid.
     */
    template <class Grid>
    void writeGrid(GridManager<Grid> const& gridManager, std::string const& gridName) const {
      writeGrid(gridManager.grid(), gridName);
    }

    /**
     * @brief writes grid into file with name gridName.
     *
     * Does not work for UGGrid.
     */
    template <class Grid>
    void writeGrid(Grid const& grid, std::string const& gridName) const {
      boost::filesystem::path path = directory/(gridName+".grid");
      boost::filesystem::path tmpPath = path;
      tmpPath += std::string("~");
      // unfortunately backup with AluGrid does not tell us if action was succesfull or not (e.g. due to wrong path)
      Dune::BackupRestoreFacility<Grid>::backup(grid, tmpPath.string());
      boost::filesystem::rename(tmpPath, path);
    }

    /**
     * @brief reads grid from file with name gridName.
     *
     * Does not work for UGGrid.
     *
     * @return pointer to the created grid or nullptr if file did not exist.
     */
    template <class Grid>
    std::unique_ptr<Grid> readGrid(std::string const& gridName) const {
      // Todo: throw if restore returns nullptr?
      return std::unique_ptr<Grid>(Dune::BackupRestoreFacility<Grid>::restore((directory/(gridName+".grid")).string()));
    }

    /**
     * @brief writes the refinement history of a grid.
     *
     * This is e.g. usefull for check pointing UGGrid, since writeGrid does not work for it.
     *
     * @param refHist is the refinement history. Suggestion for interpretation:
     * Each entry refers to one refinement step. If an entry is an empty vector
     * this means global refinement (once). Otherwise the vector for a refinement step containts pairs of refinement count
     * and cell index (similiar to the usage in the mark() method of the Dune grid interface).
     */
    template <class Index>
    void writeRefinementHistory(std::vector<std::vector<std::pair<int, Index>>> const& refHist, std::string const& gridName) const {
      boost::filesystem::path path = directory/(gridName+".refHist");
      boost::filesystem::path tmpPath = path;
      tmpPath += "~";
      std::ofstream out(tmpPath.string(), std::ios::binary);
      for(auto const& refStep : refHist) {
        out << refStep.size() << " ";
        out.write(reinterpret_cast<char const*>(refStep.data()), sizeof(std::pair<int, Index>) * refStep.size());
      }
      boost::filesystem::rename(tmpPath, path);
    }

    /**
     * @brief appends data of one refinement step of a grid to file.
     *
     * This is e.g. usefull for check pointing UGGrid, since writeGrid does not work for it.
     * The advantage compared to writeRefinementsHistory is, that if you do a lot of check-pointing and refinement steps,
     * you do not have write the whole history in each step again and again but only to append the update from the last step.
     */
    template <class Index>
    void appendRefinementHistory(std::vector<std::pair<int, Index>> const& refStep, std::string const& gridName) const {
      boost::filesystem::path path = directory/(gridName+".refHist");
      boost::filesystem::path tmpPath = path;
      tmpPath += "~";
      if(boost::filesystem::exists(path)) boost::filesystem::copy_file(path, tmpPath);
      std::ofstream out(tmpPath.string(), std::ios::binary | std::ios::app);
      out << refStep.size() << " ";
      out.write(reinterpret_cast<char const*>(refStep.data()), sizeof(std::pair<int, Index>) * refStep.size());
      boost::filesystem::rename(tmpPath, path);
    }

    /**
     * @brief reads the refinement history of a grid.
     *
     * This is e.g. usefull for check pointing UGGrid, since writeGrid does not work for it.
     * After reading in, the initial grid can be refined step by step in the same manner.
     *
     * @param refHist is the refinement history. Suggestion for interpretation:
     * Each entry refers to one refinement step. If an entry is an empty vector
     * this means global refinement (once). Otherwise the vector for a refinement step containts pairs of refinement count
     * and cell index (similiar to the usage in the mark() method of the Dune grid interface).
     */
    template <class Index>
    std::vector<std::vector<std::pair<int, Index>>> readRefinementHistory(std::string const& gridName) const {
      using RefHist = std::vector<std::vector<std::pair<int, Index>>>;
      using RefStep = typename RefHist::value_type;
      RefHist refHist;
      std::ifstream in((directory/(gridName+".refHist")).string(), std::ios::binary);
      typename RefStep::size_type nRefCells;
      while(in >> nRefCells) {
        refHist.push_back(RefStep(nRefCells));
        in.get(); // skip one white space
        in.read(reinterpret_cast<char *>(refHist.back().data()), sizeof(typename RefStep::value_type) * nRefCells);
      }
      return refHist;
    }

    /**
     * @brief writes data into file with name dataName. Version for BlockVector.
     *
     * @tparam Array is template type, which behaves like an array, e.g. Dune::FieldVector
     */
    template <template<class, int> class Array, class Scalar, int dim>
    void writeData(Dune::BlockVector<Array<Scalar, dim>> const& data, std::string const& dataName) const {
      boost::filesystem::path path = directory/(dataName+".data");
      boost::filesystem::path tmpPath = path;
      tmpPath += std::string("~");
      writeData(data, std::ofstream(tmpPath.string(), std::ios::binary));
      boost::filesystem::rename(tmpPath, path);
    }

    /**
     * @brief writes data into file with name dataName. Version for LinearProductSpace (of FE-functions or coefficient vectors).
     *
     * @tparam DataBlocks has to be something like LinearProductSpace, i.e. e.g. VariableSet or CoefficientVectorRepresentation.
     */
    template <class DataBlocks>
    void writeData(DataBlocks const& data, std::string const& dataName) const {
      boost::filesystem::path path = directory/(dataName+".data");
      boost::filesystem::path tmpPath = path;
      tmpPath += std::string("~");
      std::ofstream out(tmpPath.string(), std::ios::binary);
      boost::fusion::for_each(data.data, [this, &out](auto const& dataBlock) {
        out = this->writeData(dataBlock, std::move(out));
      });
      boost::filesystem::rename(tmpPath, path);
    }

    /**
     * @brief writes data into file with name dataName. Version for FE-function.
     */
    template <class FeSpace, int m>
    void writeData(FunctionSpaceElement<FeSpace, m> const& data, std::string const& dataName) const {
      writeData(data.coefficients(), dataName);
    }

    /**
     * @brief reads data from file with name dataName. Version for BlockVector.
     *
     * @tparam Array is template type, which behaves like an array, e.g. Dune::FieldVector
     */
    template <template<class, int> class Array, class Scalar, int dim>
    void readData(Dune::BlockVector<Array<Scalar, dim>> & data, std::string const& dataName) const {
      readData(data, std::ifstream((directory/(dataName+".data")).string(), std::ios::binary));
    }

    /**
     * @brief reads data from file with name dataName. Version for LinearProductSpace (of FE-functions or coefficient vectors).
     *
     * @tparam DataBlocks has to be something like LinearProductSpace, i.e. e.g. VariableSet or CoefficientVectorRepresentation.
     */
    template <class DataBlocks>
    void readData(DataBlocks & data, std::string const& dataName) const {
      std::ifstream in((directory/(dataName+".data")).string(), std::ios::binary);
      boost::fusion::for_each(data.data, [this, &in](auto & dataBlock) {
        in = this->readData(dataBlock, std::move(in));
      });
    }

    /**
     * @brief reads data from file with name dataName. Version for FE-function.
     */
    template <class FeSpace, int m>
    void readData(FunctionSpaceElement<FeSpace, m> & data, std::string const& dataName) const {
      readData(data.coefficients(), dataName);
    }

  protected:
    void setTimeStamp();

    template <template<class, int> class Array, class Scalar, int dim>
    std::ofstream writeData(Dune::BlockVector<Array<Scalar, dim>> const& data, std::ofstream out) const {
      out.write(reinterpret_cast<char const*>(&data[0]), sizeof(Scalar) * data.size() * dim);
      return out;
    }

    template <class FeSpace, int m>
    std::ofstream writeData(FunctionSpaceElement<FeSpace, m> const& data, std::ofstream out) const {
      return writeData(data.coefficients(), std::move(out));
    }

    template <template<class, int> class Array, class Scalar, int dim>
    std::ifstream readData(Dune::BlockVector<Array<Scalar, dim>> & data, std::ifstream in) const {
      in.read(reinterpret_cast<char *>(&data[0]), sizeof(Scalar) * data.size() * dim);
      return in;
    }

    template <class FeSpace, int m>
    std::ifstream readData(FunctionSpaceElement<FeSpace, m> & data, std::ifstream in) const {
      return readData(data.coefficients(), std::move(in));
    }

    std::string timeStamp_;
    boost::filesystem::path directory;
  };
}
#endif // CHECK_POINT_HH
