#pragma once
#include <iomanip>
#include <chrono>
#include <vector>
#include <string>



// JFG 2022: This is a direct copy of a bunch of stuff from my JSL library which I have copied verbatim here so that the entirety of JSL doesn't need to be copied across. 

namespace MCMC
{
	/*!
		Given a duration in seconds, convert it into standard Day/Hour/Minute/Second formatted string. Times less than 1 second are reported as "less than 1 second"
		\param timeInSeconds The time to be converted
		\returns A human-readable string equal to the input
	 */
	inline std::string FormatDuration(int timeInSeconds)
	{
		
		std::vector<std::string> divisions= {"Day", "Hour", "Minute", "Second"};
		std::vector<int> duration = {86400, 3600, 60 , 1};
		//std::vector<int> chunks = std::vector(0,divisions.size());
		
		bool nothingAdded = true;
		std::string output = "";
		if (timeInSeconds < 0)
		{
			timeInSeconds = abs(timeInSeconds);
			output += "(negative) ";
		}
		for (int i = 0; i < divisions.size(); ++i)
		{
			int v = timeInSeconds / duration[i];
			timeInSeconds -= v * duration[i];
			
			if (v > 0)
			{
				nothingAdded = false;
				output += std::to_string(v) + " " + divisions[i];
				if (v > 1)
				{
					output += "s";
				}
				output += " ";
			}
		}
		if (nothingAdded)
		{
			output = "less than 1 second";
		}	
		return output;
	}
	



	//!Deletes the last linebreak character, jumping the cursor up one line
	inline void jumpLineUp()
	{
		printf("\033[F");
	}
	//!Equivalent to calling the "clear" command on the shell, removes all text on the terminal, and returns the cursor to the home position.
	inline void clearScreen()
	{
		printf("\033[H\033[J\033[F");
	}
	//!Removes the last line of text, moving the cursor from its current position, to the beginning of the line. If you want to remove the linebreak which caused that line to exist in the first place, must be followed up with jumpLineUp(). ALternatively, if a linebreak has been printed and you want to delete that line *first*, you must call jumpLineUp() before deleteLine().
	inline void deleteLine()
	{
		printf("\33[2K\r");
		// std::cout << std::flush;
	}
	/*!
		A neat display widget for tracking the progress of a process by printing a progress bar on the screen, which updates to visually indicate the process is running. Supports arbitrary "levels" of progress bars which can be independently modified through a variadic argument interface.
		\tparam Dimension The number of concurrent bars to print (fixes the number of arguments which must be provided to subsequent functions). Defaults to 1 if 0 template arguments provided
		\tparam DeleteMode If true, utilises JSL::deleteLine() to remove the previous line before reprinting, such that the bar appears animated. If false, provides a more limited form of animation (false only available in Dimension==1 mode). Defaults to true.
		\tparam Symbol The symbol which is printed within the bar to fill up space. Defaults to "#".
		\tparam MaxHashes The number of symbols in a full bar (dictates the width of the progress bar)
	*/
	template<int Dimension = 1, bool DeleteMode = true, char Symbol = '#',unsigned int MaxHashes =20>
	class ProgressBar
	{
		protected:
			bool firstPrint = true; //!<Toggle to ensure open-brace "[" is printed only once
			std::vector<int> BarTargets; //!<Storage for the target number of processes
			std::vector<double> BarProgress; //!< Current progress value (between 0 and 1)
			std::vector<int> Hashes; //!< Current printed number of symbols for each bar
			std::vector<std::string> Names; //!< The names printed in front of each bar 

			bool reprintNeeded; //!<Toggle which ensures reprinting only occurs when something changes, prevents unnecessary overhead
			int bufferWidth = 10; //!<Default size allocated for names, increases if a value of anmes would overfill it
			int prevHashes = 0; //!<Used in DeleteMode=false to compute new hashe number

			//! The variadic internal function called by the compiler, looping over the arbitrary number of inputs and storing them in internal vectors
			template<typename... Ts>
			void UnpackTargets(int target, Ts... targets)
			{
				
				BarTargets.push_back(target);
				BarProgress.push_back(0);
				Hashes.push_back(0);
				if constexpr (sizeof...(targets) > 0)
				{
					UnpackTargets(targets...);
				}
			}

			//! The variadic internal function which is called by Update(), loops over the arguments and computes the number of Symbols required by each bar
			template<typename... Ts>
			void UnrollPositions(const int idx,int pos, Ts... remainder)
			{
				BarProgress[idx]  = (double)pos/(BarTargets[idx]-1);
				int newHashes = std::min(BarProgress[idx],1.0) * MaxHashes;
				
				if (newHashes != Hashes[idx])
				{
					prevHashes = Hashes[idx];
					reprintNeeded = true;
					Hashes[idx] = newHashes;
				}
				if constexpr( sizeof...(remainder)>0)
				{
					UnrollPositions(idx+1,remainder...);
				}
				
				
			}

			//! Loops over the existing bars and erases them using ANSI codes
			void ClearScreen()
			{
				for (int i = 0; i < Dimension; ++i)
				{
					jumpLineUp();
					deleteLine();
				}
			}
			//! The version of Printing which happens if DeleteMode is true -- clears the previous bars using ANSI codes, then overwrites them with the new values
			virtual void PrintPositions_DeleteMode()
			{
				if (!firstPrint)
				{
					ClearScreen();
				}
				firstPrint = false;
				for (int i = 0; i < Dimension; ++i)
				{
					int newHashes = Hashes[i];
					if (Names.size() > 0)
					{
						std::cout << std::setw(bufferWidth)  << Names[i]<< " ";
					}
					std::cout << "[" << std::string(newHashes,Symbol) << std::string(MaxHashes - newHashes,' ') << "]\n";
				}
			}

			//! The version of Printing which happens if DeleteMode is false - gradually prints the progress bar one Symbol at a time, then adds the close brace "]" whn the final one is printed
			void PrintPositions_RetainMode()
			{
		
				int newHashes = Hashes[0] - prevHashes;
				if (prevHashes == 0)
				{
					if (Names.size() == 1)
					{
						std::cout << std::setw(bufferWidth)  << Names[0]<< " ";
					}
					std::cout << "[";
				}
				std::cout << std::string(newHashes,Symbol);
				if (Hashes[0] == MaxHashes)
				{
					std::cout << "]";
				}
				std::cout << std::flush;
			}
		public:

			/*!
				Constructor function. Accepts a number of integer arguments equal to Dimension, throws an error if this is not true.
				\param targetCounts The number of processes (i.e. max loop index) that must complete before each bar is considered full.
			*/
			template<typename... Ts>
			ProgressBar(Ts... targetCounts) 
			{
				static_assert(sizeof...(targetCounts) == Dimension);
				static_assert(Dimension > 0);
				static_assert(MaxHashes > 1);
				UnpackTargets(targetCounts...);
			}

			/*!
				Provide a progress update (i.e. loop index) for each of the active bars. Accepts a number of integer arguments equal to Dimension.
				\param positons The current process count (to be compared against the targets) of the open progress bars
			*/
			template<typename... Ts>
			void Update(Ts... positions)
			{
				static_assert(sizeof...(positions) == Dimension);
				reprintNeeded  = false;
				UnrollPositions(0,positions...);

				//only reprint the bar if something has changed
				if (reprintNeeded)
				{
					if constexpr (DeleteMode)
					{
						PrintPositions_DeleteMode();
					}
					else
					{
						PrintPositions_RetainMode();
					}
				}
			}

			

			//! Sets all names for the bars simultaneously \param names A vector of strings (equal to Dimension), each of which is printed before their respective bar.
			void SetName(const std::vector<std::string> & names)
			{
				assert(names.size() == Dimension);//"Name vector must be same size as bar count"
				Names = names;
				for (int i = 0; i < Dimension; ++i)
				{
					bufferWidth = std::max(bufferWidth,(int)names[i].length());
				}
			}

			//! Changes the name of the idx-th bar \param idx The bar name to be changes \param name The string which is printed in the front of the bar
			void SetName(unsigned int idx, const std::string & name)
			{
				assert(idx < Dimension);//"Assign name to bar between 0 and BarCount-1"
				if (Names.size() < Dimension)
				{
					Names = std::vector<std::string>(Dimension,"");
				}
				Names[idx] = name;
				bufferWidth = std::max(bufferWidth,(int)name.length());
			}
			
			//!An override for SetName(unsigned int idx, const std::string & name) available only when Dimension == 1 \param name The message to go in front of the bar as it prints
			void SetName(const std::string & name)
			{
				static_assert(Dimension == 1);
				Names = {name};
				bufferWidth = std::max(bufferWidth,(int)name.length());
			}

			//! Deletes all levels of the bar from the screen, useful once the bar has finished and is no longer needed
			void Clear()
			{
				ClearScreen();
			}
	};


	template<int Dimension = 1, char Symbol = '#',unsigned int MaxHashes =20>
	class PredictionBar : public ProgressBar<Dimension,true,Symbol,MaxHashes>
	{
		public:

			template<typename... Ts>
			PredictionBar(Ts... targetCounts) : ProgressBar<Dimension,true,Symbol,MaxHashes>(targetCounts...)
			{
				Start = std::chrono::system_clock::now();
				UpdateTime = Start;
				// Update(0);
			}


			template<typename... Ts>
			void Update(Ts... positions)
			{
				static_assert(sizeof...(positions) == Dimension);
				this->reprintNeeded  = false;
				this->UnrollPositions(0,positions...);
				ComputeTime();
				//only reprint the bar if something has changed
				if (this->reprintNeeded)
				{
					PrintPositions_DeleteMode();

				}
			}

		private:
			std::chrono::time_point<std::chrono::system_clock> Start;
			std::chrono::time_point<std::chrono::system_clock> UpdateTime;
			double CurrentDuration;
			double Prediction = 0;
			double prevProgress = 0;
			void ComputeTime()
			{
				std::chrono::time_point<std::chrono::system_clock> Now = std::chrono::system_clock::now();
				std::chrono::duration<double,std::ratio<1,1>> totalDuration = Now - Start;
				std::chrono::duration<double,std::ratio<1,1>> updateDuration = Now - UpdateTime;
				CurrentDuration = totalDuration.count();
				if (updateDuration.count() > 1)
				{
					this->reprintNeeded = true;
				}
				if (this->reprintNeeded)
				{
					UpdateTime = Now; // also catches reprints due to altered barstates!
				
				
					double totalProg = 0;
					double next = 1;
					for (int i = 0; i < Dimension; ++i)
					{
						double prog = this->BarProgress[i] * (this->BarTargets[i]-1)/(this->BarTargets[i]);
						totalProg += prog * next;
						next = next/(this->BarTargets[i]);
						// std::cout << i << "  " << prog << "  " << next << "  " << totalProg << std::endl;
					}
					
					double globalRate = totalProg/totalDuration.count();
					double stepRate = (totalProg - prevProgress)/updateDuration.count();
					double w = 0.2;
					double meanRate = (1.0 - w)*globalRate + w*stepRate;
					double anticipate = std::max(0.0,Prediction - updateDuration.count());
					double mem = 0.9;
					double activePrediction = (1.0 - totalProg)/meanRate;
					if (activePrediction > anticipate)
					{
						mem = 0.5;
					}
					
					// std::cout << globalRate << "  " << stepRate << std::endl;
					// std::cout << "Anticiapting " << Prediction << " - " << updateDuration.count() << " = " << Prediction - updateDuration.count() << " = " << anticipate << " actively " << activePrediction << " av = " << mem * anticipate + (1.0 - mem)*activePrediction<< "\n\n\n\n";
					Prediction = mem * anticipate + (1.0 - mem)*activePrediction;
					// std::cout << "Estimated " << totalProg << " through at rates " << globalRate << "   " << stepRate << "\n\n\n\n";
					prevProgress = totalProg;
				}
			}
			virtual void PrintPositions_DeleteMode()
			{

				if (!this->firstPrint)
				{
					this->ClearScreen();
				}
				this->firstPrint = false;
				for (int i = 0; i < Dimension; ++i)
				{
					int newHashes = this->Hashes[i];
					if (this->Names.size() > 0)
					{
						std::cout << std::setw(this->bufferWidth)  << this->Names[i]<< " ";
					}
					std::cout << "[" << std::string(newHashes,Symbol) << std::string(MaxHashes - newHashes,' ') << "]";
					if (i == Dimension - 1)
					{
						std::cout << " ETR: " << FormatDuration(Prediction);
					}
					// if (i == 0)
					// {
					// 	std::cout << " " << (int)(100*prevProgress) << "%" << " in " << JSL::FormatDuration(CurrentDuration);
					// }
					std::cout << "\n";
				}
			}
	};


}