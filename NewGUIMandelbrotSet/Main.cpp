#include <SFGUI/SFGUI.hpp>
#include <SFGUI/Widgets.hpp>
#include <SFML/Graphics.hpp>
#include <string>
#include <iostream>
#include <thread>
#include <ctime>
#include <mpir.h>
#include <mpirxx.h>
#include <vector>
#include <complex>
#include <cmath>
#include <Windows.h>

enum colortype { GRAYSCALE, HUECYCLE };

class Set
{
private:
	bool IsGenerated = false;
	int completedThreads = 0;

	std::complex<double>* pixels;
	int* pixelsIter;

	std::vector<std::thread> threads;
	std::vector<std::complex<double>> orbit;

	void CalculateOrbit()
	{
		mpf_class x(0,1000);
		mpf_class y(0,1000);
		mpf_class xtemp;

		for (int i = 0; i < maxIter; i++)
		{
			orbit.push_back(std::complex<double>(x.get_d(), y.get_d())*2.0); //Предумножение на 2 для оптимизации
			xtemp = x*x - y*y + orbitPos.real();
			y = 2.0 * x*y + orbitPos.imag();
			x = xtemp;

			if (x > 1024 || x < -1024 || y > 1024 || y < -1024)
			{
				return;
			}
		}
	}

	void ThreadTask(int threadNumber, int threadsNumber)
	{
		mpf_class x0(0,1000), y0(0,1000);
		double halfx = (double)width / 2;
		double halfy = (double)height / 2;
		for (double i = threadNumber; i < width; i += threadsNumber) //Потоки рассчитывают пиксели по отдельным столбам. 0 столб - 0 поток, 1 столб - 1 поток и тд.
		{
			for (double j = 0; j < height; j++)
			{
				x0 = (i - halfx) / width * zoom + xpos;
				y0 = -(j - halfy) / width * zoom + ypos;

				//pixels[height*(int)i + (int)j] = CalculateSSAAPixel(x0, y0);
				pixelsIter[height*(int)i + (int)j] = CalculateSSAAPixel(x0, y0);
			}
		}

		completedThreads++;
	}

	//std::complex<double> CalculateSSAAPixel(mpf_class& x, mpf_class& y)
	int CalculateSSAAPixel(mpf_class& x, mpf_class& y)
	{
		int* subpixelsIter = new int[SSAALevelSquared];
		int orbitIterNumber = (int)orbit.size();

		mpf_class tempx(0,1000), tempy(0,1000);

		for (int i = 0; i < SSAALevel; i++) //Calculate all subpixels
		{
			for (int j = 0; j < SSAALevel; j++)
			{
				tempx = x + SSAAStep*i - orbitPos.real();
				tempy = y + SSAAStep*j - orbitPos.imag();
				std::complex<double> d0(tempx.get_d(), tempy.get_d());
				std::complex<double> dn = d0;

				int iter = 1;

				while (iter < orbitIterNumber && std::norm(dn + orbit[iter]*0.5) < 4)
				{
					dn = orbit[iter] * dn + dn*dn + d0; //Умножение на 2 перенесено в CalculateOrbit
					iter++;
				}

				int pos = SSAALevel*i + j;
				subpixelsIter[pos] = iter;
			}
		}

		int avgIter = 0;

		for (int i = 0; i < SSAALevelSquared; i++)
		{
			avgIter += subpixelsIter[i];
		}
		avgIter /= SSAALevelSquared;

		delete[] subpixelsIter;

		return avgIter;
	}

	sf::Color HSVToRGB(double hue, double saturation, double brightness) //0-360, 0-255, 0-255
	{
		int hi = (int)(hue / 60) % 6;
		int Vmin = (int)std::round(((255 - saturation) * brightness) / 255);
		double a = ((brightness - Vmin) * ((std::remainder(hue, 60.0)) / 60));
		int Vinc = (int)std::round((Vmin + a));
		int Vdec = (int)std::round(brightness - a);
		switch (hi)
		{
		case 0:
			return sf::Color((int)brightness, Vinc, Vmin);
		case 1:
			return sf::Color(Vdec, (int)brightness, Vmin);
		case 2:
			return sf::Color(Vmin, (int)brightness, Vinc);
		case 3:
			return sf::Color(Vmin, Vdec, (int)brightness);
		case 4:
			return sf::Color(Vinc, Vmin, (int)brightness);
		case 5:
			return sf::Color((int)brightness, Vmin, Vdec);
		}
		return sf::Color((int)brightness, Vmin, Vdec);
	}

	mpf_class xpos;
	mpf_class ypos;
	int width;
	int height;
	int maxIter;
	int SSAALevel;
	int SSAALevelSquared;
	std::complex<mpf_class> orbitPos;

	mpf_class zoom;
	mpf_class SSAAStep;

public:


	Set(mpf_class& xp, mpf_class& yp, mpf_class& z, mpf_class& orbx, mpf_class& orby, int mI, int ssaa, int w, int h) : xpos(xp), ypos(yp), maxIter(mI), SSAALevel(ssaa), width(w), height(h), zoom(z) {
		xpos.set_prec(1000);
		ypos.set_prec(1000);
		zoom.set_prec(1000);
		SSAAStep.set_prec(1000);

		pixelsIter = new int[width*height];
		SSAALevelSquared = SSAALevel*SSAALevel;
		SSAAStep = (zoom / width) / SSAALevel;

		orbitPos = std::complex<mpf_class>(orbx, orby);
	}

	void Generate(int threadsNumber)
	{
		if (IsGenerated)
		{
			return;
		}
		IsGenerated = true;

		CalculateOrbit();
		for (int i = 0; i < threadsNumber; i++)
		{
			threads.push_back(std::thread(&Set::ThreadTask, this, i, threadsNumber));
		}
		for (int i = 0; i < threadsNumber; i++)
		{
			threads[i].detach();
		}
	}

	bool IsGenerating()
	{
		return !(completedThreads == threads.size());
	}

	sf::Image GetImage(colortype type)
	{
		if (IsGenerating() || !IsGenerated)
		{
			return sf::Image();
		}

		sf::Image image;
		image.create(width, height, sf::Color::Black);

		if (type == HUECYCLE)
		{
			for (int i = 0; i < width; i++)
			{
				for (int j = 0; j < height; j++)
				{
					//double nu = iter - std::log2(std::log2(zn_size));
					//int i = (int)(nu * 10) % gradient.size();

					//return gradient[i];

					//double clr = std::norm(pixels[height*i + j]);
					int iter = pixelsIter[height*i + j];
					image.setPixel(i, j, HSVToRGB(iter, 255, 255));
				}
			}
		}
		else if (type == GRAYSCALE)
		{
			for (int i = 0; i < width; i++)
			{
				for (int j = 0; j < height; j++)
				{
					int iter = pixelsIter[height*i + j];
					image.setPixel(i, j, sf::Color(iter, iter, iter));
				}
			}
		}
		return image;
	}

	~Set()
	{
		delete[] pixels;
		delete[] pixelsIter;
	}
};

class GUI
{
private:
	std::string dtos(double d)
	{
		std::stringstream str;
		str << std::fixed << std::setprecision(15) << d;
		std::string strn = str.str();
		return strn;
	}

	std::string mpftostr(mpf_class& d)
	{
		std::stringstream str;
		str.setf(std::stringstream::scientific);
		str.precision(100);
		str << d;
		std::string strn = str.str();
		return strn;
	}

	void controlThreads()
	{
		if (curSet != nullptr)
		{
			if (!curSet->IsGenerating() && calculating)
			{
				std::string time = std::to_string(std::round(clock.getElapsedTime().asSeconds() * 100) / 100);
				time.erase(time.begin() + time.find(".") + 3, time.end());

				generateTimeLabel->SetText("Finished in " + time + " seconds");
				calculating = false;
				spinner->Stop();

				if (autoMakeCheckButton->IsActive())
				{
					exportbut();
				}
			}
		}
	}

	/*
	void generatePart(int part, int partNum, int width, int height, int* arr)
	{
		double xInput = std::stod(xEntry->GetText().toAnsiString());
		double yInput = std::stod(yEntry->GetText().toAnsiString());
		double zoomInput = std::stod(zoomEntry->GetText().toAnsiString());

		int maxIter = std::stoi(maxIterEntry->GetText().toAnsiString());
		//int width = std::stoi(widthEntry->GetText().toAnsiString());
		//int height = std::stoi(heightEntry->GetText().toAnsiString());

		double halfx = (double)width / 2;
		double halfy = (double)height / 2;

		double x0, y0, x, y, xtemp;
		int iteration, pos;

		double zoomInputDivHalfX = zoomInput / halfx;
		double minusZoomInputDivHalfX = -zoomInputDivHalfX;
		double minusZoomInputPlusXInput = -zoomInput + xInput;
		double zoomInputDivHalfXTimesHalfY = zoomInputDivHalfX * halfy;


		for (double i = part; i < width; i += partNum)
		{
			for (double j = 0; j < height; j++)
			{
				//x0 = (i - halfx) / halfx * zoomInput + xInput; //Было
				//y0 = -(j - halfy) / halfx * zoomInput + yInput;

				x0 = i * zoomInputDivHalfX + minusZoomInputPlusXInput; //Стало
				y0 = j * minusZoomInputDivHalfX + zoomInputDivHalfXTimesHalfY + yInput;


				x = 0.0;
				y = 0.0;

				iteration = 0;

				double xsq = 0.0;
				double ysq = 0.0;

				while (xsq + ysq < 4 && iteration < maxIter) {
					xtemp = xsq - ysq + x0;
					y = 2 * x*y + y0;
					x = xtemp;
					iteration++;

					xsq = x*x;
					ysq = y*y;
				}

				pos = (int)i * height + (int)j;

				arr[pos] = iteration;
				//arr[pos] = (int)(x/y*iteration);
			}
		}
		finishedThreads++;
	}

	void generatePartSuperSampled(int part, int partNum, int width, int height, int supersampling, int* arr)
	{
		double xInput = std::stod(xEntry->GetText().toAnsiString());
		double yInput = std::stod(yEntry->GetText().toAnsiString());
		double zoomInput = std::stod(zoomEntry->GetText().toAnsiString());

		int maxIter = std::stoi(maxIterEntry->GetText().toAnsiString());
		//int width = std::stoi(widthEntry->GetText().toAnsiString());
		//int height = std::stoi(heightEntry->GetText().toAnsiString());

		double halfx = (double)width / 2;
		double halfy = (double)height / 2;

		double sampStep = 1 / (double)supersampling;

		int sampleNum = supersampling*supersampling;
		int* numbers = new int[sampleNum]; //части одного пикселя

		for (double i = part; i < width; i += partNum)
		{
			for (double j = 0; j < height; j++)
			{

				int num = 0; //номер части пикселя в массиве

				for (double si = 0; si < 1; si += sampStep)
				{
					for (double sj = 0; sj < 1; sj += sampStep)
					{

						double x0 = (i + si - halfx) / halfx * zoomInput + xInput;
						double y0 = -(j + sj - halfy) / halfx * zoomInput + yInput;

						double x = 0.0;
						double y = 0.0;

						int iteration = 0;

						while (x*x + y*y < 4 && iteration < maxIter) {
							double xtemp = x*x - y*y + x0;
							y = 2 * x*y + y0;
							x = xtemp;
							iteration++;
						}
						numbers[num] = iteration;
						num++;

					}
				}

				int pos = (int)i * height + (int)j;
				int avg = 0;
				for (int i = 0; i < sampleNum; i++)
				{
					avg += numbers[i];
				}
				avg /= sampleNum;

				arr[pos] = avg;
			}
		}
		delete[] numbers;

		finishedThreads++;
	}

	*/

	void generate()
	{
		spinner->Start();
		clock.restart();

		calculating = true;

		GetInputs();

		if (curSet != nullptr)
		{
			delete curSet;
		}
		curSet = new Set(inputX, inputY, inputZoom, inputOrbitX, inputOrbitY, inputIter, inputSSAA, inputWidth, inputHeight);
		curSet->Generate(inputThreadsNum);

	}

	void estimateTime()
	{
		sf::Clock timeclock;
		timeclock.restart();

		GetInputs();

		int width = 256;
		int height = 256;

		Set set(inputX, inputY, inputZoom, inputOrbitX, inputOrbitY, inputIter, inputSSAA, width, height);

		while (set.IsGenerating())
		{
			std::cout << "zz";
			//BAD CODE
		}

		double time = (double)timeclock.getElapsedTime().asSeconds();

		time *= (inputWidth * inputHeight) / (width*height);
		std::string timeS = std::to_string(std::round(time * 100) / 100);
		timeS.erase(timeS.begin() + timeS.find(".") + 3, timeS.end());
		estimatedTimeLabel->SetText("Est. time: " + timeS);
	}

	void exportbut()
	{
		if (curSet == nullptr)
		{
			return;
		}

		GetInputs();

		sf::Image image;
		image = curSet->GetImage(inputType);
		image.saveToFile(filenameEntry->GetText().toAnsiString() + ".png");

		if (autoMakeCheckButton->IsActive())
		{
			inputZoom *= inputZoomStep;
			SetInputs();
			std::string text = filenameEntry->GetText().toAnsiString();
			filenameEntry->SetText(std::to_string(std::time(0)));
			generate();
		}
	}

	std::string parseDoubleString(std::string x)
	{
		for (int i = (int)x.size() - 1; i >= 0; i--)
		{
			char ch = x[i];
			if ((ch < '0' || ch > '9') && ch != '.' && ch != ',' && ch != '-' && ch != 'e')
			{
				x.erase(i, 1);
			}
			if (ch == ',')
			{
				x[i] = '.';
			}
		}
		if (x.size() == 0)
		{
			x.push_back('0');
		}
		return x;
	}

	std::string parsePositiveDoubleString(std::string x)
	{
		for (int i = (int)x.size() - 1; i >= 0; i--)
		{
			char ch = x[i];
			if (ch == 45)
			{
				x.erase(i, 1);
			}
		}
		return parseDoubleString(x);
	}

	std::string parseIntString(std::string x)
	{
		for (int i = (int)x.size() - 1; i >= 0; i--)
		{
			char ch = x[i];
			if ((ch < 47 || ch > 57))
			{
				x.erase(x.begin() + i);
			}
		}
		if (x.size() == 0)
		{
			x.push_back('0');
		}
		return x;
	}

	void parseInputs()
	{
		xEntry->SetText(parseDoubleString(xEntry->GetText().toAnsiString()));
		yEntry->SetText(parseDoubleString(yEntry->GetText().toAnsiString()));
		zoomEntry->SetText(parseDoubleString(zoomEntry->GetText().toAnsiString()));
		widthEntry->SetText(parseIntString(widthEntry->GetText().toAnsiString()));
		heightEntry->SetText(parseIntString(heightEntry->GetText().toAnsiString()));
		maxIterEntry->SetText(parseIntString(maxIterEntry->GetText().toAnsiString()));
		zoomStepEntry->SetText(parseDoubleString(zoomStepEntry->GetText().toAnsiString()));
	}

	void calcPreview() //Рассчитывает и выводит превью
	{
		GetInputs();

		Set preview(inputX, inputY, inputZoom, inputOrbitX, inputOrbitY, inputIter, 1, 400, 400);
		preview.Generate(inputThreadsNum);

		while (preview.IsGenerating())
		{
			std::cout << "gen";
			//REDO, BAD CODE!
		}

		auto pic = sfg::Image::Create();
		pic->SetImage(preview.GetImage(inputType));

		std::vector<std::shared_ptr<sfg::Widget>> children = previewBox->GetChildren();
		if (children.size() > 1)
		{
			previewBox->Remove(children[1]);
		}
		previewBox->PackStart(pic);
	}

	void GetInputs()
	{
		parseInputs();

		inputX.set_str(xEntry->GetText().toAnsiString(), 10);
		inputY.set_str(yEntry->GetText().toAnsiString(), 10);
		inputZoom.set_str(zoomEntry->GetText().toAnsiString(), 10);

		switch (SSAAComboBox->GetSelectedItem())
		{
		case 0:
			inputSSAA = 1;
			break;
		case 1:
			inputSSAA = 2;
			break;
		case 2:
			inputSSAA = 4;
			break;
		case 3:
			inputSSAA = 8;
			break;
		case 4:
			inputSSAA = 16;
			break;
		}

		inputIter = std::stoi(maxIterEntry->GetText().toAnsiString());
		inputWidth = std::stoi(widthEntry->GetText().toAnsiString());
		inputHeight = std::stoi(heightEntry->GetText().toAnsiString());
		inputZoomStep = std::stoi(zoomStepEntry->GetText().toAnsiString());

		if (colorComboBox->GetSelectedItem() == 1) //Hue cycle
		{
			inputType = HUECYCLE;
		}
		else
		{
			inputType = GRAYSCALE;
		}

		inputThreadsNum = (int)threadNumSpinbutton->GetValue();
	}

	void SetInputs()
	{
		//mp_exp_t exp;
		//xEntry->SetText(inputX.get_str(exp));
		//yEntry->SetText(inputY.get_str(exp));
		//zoomEntry->SetText(inputZoom.get_str(exp));
		xEntry->SetText(mpftostr(inputX));
		yEntry->SetText(mpftostr(inputY));
		zoomEntry->SetText(mpftostr(inputZoom));
	}

	mpf_class inputX;
	mpf_class inputY;
	mpf_class inputZoom;
	int inputIter;
	int inputSSAA;
	int inputWidth;
	int inputHeight;
	int inputThreadsNum;
	colortype inputType;
	double inputZoomStep;

	mpf_class inputOrbitX = 0;
	mpf_class inputOrbitY = 0;

	Set* curSet;

public:
	int *setArray;
	int genwidth = 0, genheight = 0;

	int finishedThreads = 0;
	int threadsAtStart = 0;
	bool calculating = false;
	sf::Clock clock; //Для расчета времени генерации

	sf::RenderWindow app_window;
	sfg::SFGUI sfgui;
	sfg::Desktop desktop;

	sfg::Window::Ptr settingsWindow;
	sfg::Window::Ptr previewWindow;
	sfg::Window::Ptr exportWindow;

	sfg::CheckButton::Ptr autoMakeCheckButton;

	sfg::Label::Ptr xLabel;
	sfg::Label::Ptr yLabel;
	sfg::Label::Ptr maxIterLabel;
	sfg::Label::Ptr zoomLabel;
	sfg::Label::Ptr widthLabel;
	sfg::Label::Ptr heightLabel;
	sfg::Label::Ptr colorLabel;
	sfg::Label::Ptr threadNumLabel;
	sfg::Label::Ptr generateTimeLabel;
	sfg::Label::Ptr estimatedTimeLabel;
	sfg::Label::Ptr controlsLabel;
	sfg::Label::Ptr SSAALabel;
	sfg::Label::Ptr zoomStepLabel;


	sfg::Entry::Ptr xEntry;
	sfg::Entry::Ptr yEntry;
	sfg::Entry::Ptr zoomEntry;
	sfg::Entry::Ptr maxIterEntry;
	sfg::Entry::Ptr widthEntry;
	sfg::Entry::Ptr heightEntry;
	sfg::Entry::Ptr filenameEntry;
	sfg::Entry::Ptr zoomStepEntry;

	sfg::SpinButton::Ptr threadNumSpinbutton;

	sfg::Button::Ptr previewButton;
	sfg::Button::Ptr generateButton;
	sfg::Button::Ptr exportButton;
	sfg::Button::Ptr estimateTimeButton;

	sfg::Box::Ptr settingsBox;
	sfg::Box::Ptr xBox;
	sfg::Box::Ptr yBox;
	sfg::Box::Ptr zoomBox;
	sfg::Box::Ptr maxIterBox;
	sfg::Box::Ptr previewBox;
	sfg::Box::Ptr widthBox;
	sfg::Box::Ptr heigthBox;
	sfg::Box::Ptr exportBox;
	sfg::Box::Ptr threadBox;
	sfg::Box::Ptr SSAABox;

	sfg::ComboBox::Ptr colorComboBox;

	sfg::Spinner::Ptr spinner;

	sfg::ComboBox::Ptr SSAAComboBox;

	GUI()
	{
		inputX.set_prec(100);
		inputY.set_prec(100);
		inputZoom.set_prec(100);
		inputOrbitX.set_prec(100);
		inputOrbitY.set_prec(100);
	}

	void SetupGUI()
	{
		app_window.create(sf::VideoMode(1280, 720), "Mandelbrot Set Generator", sf::Style::Titlebar | sf::Style::Close);
		// We have to do this because we don't use SFML to draw.
		app_window.resetGLStates();
		app_window.setFramerateLimit(60);


		//Настройка окон
		settingsWindow = sfg::Window::Create();
		settingsWindow->SetTitle("Settings");
		settingsWindow->SetRequisition(sf::Vector2f(1280.f, 200.f));

		previewWindow = sfg::Window::Create();
		previewWindow->SetTitle("Preview");
		previewWindow->SetPosition(sf::Vector2f(0.f, 200.f));
		previewWindow->SetRequisition(sf::Vector2f(422.f, 446.f));

		exportWindow = sfg::Window::Create();
		exportWindow->SetTitle("Export");
		exportWindow->SetPosition(sf::Vector2f(1022.f, 200.f));
		//exportWindow->SetRequisition(sf::Vector2f(200.f, 300.f));

		//Интерфейс окна настроек
		xLabel = sfg::Label::Create("x:"); //x
		xEntry = sfg::Entry::Create("-1");
		yLabel = sfg::Label::Create("y:"); //y
		yEntry = sfg::Entry::Create("0");
		zoomLabel = sfg::Label::Create("Zoom:"); //zoom
		zoomEntry = sfg::Entry::Create("2");
		maxIterLabel = sfg::Label::Create("Max iter:"); //iter
		maxIterEntry = sfg::Entry::Create("255");
		previewButton = sfg::Button::Create();
		previewButton->SetLabel("Preview");
		previewButton->GetSignal(sfg::Widget::OnLeftClick).Connect([this] {calcPreview(); });


		xEntry->SetRequisition(sf::Vector2f(1000.f, 10.f)); //Размеры текстовых полей
		yEntry->SetRequisition(sf::Vector2f(1000.f, 10.f));
		zoomEntry->SetRequisition(sf::Vector2f(1000.f, 10.f));
		maxIterEntry->SetRequisition(sf::Vector2f(1000.f, 10.f));

		settingsBox = sfg::Box::Create(sfg::Box::Orientation::VERTICAL); //Коробки для правильного расположения всего
		xBox = sfg::Box::Create(sfg::Box::Orientation::HORIZONTAL);
		yBox = sfg::Box::Create(sfg::Box::Orientation::HORIZONTAL);
		zoomBox = sfg::Box::Create(sfg::Box::Orientation::HORIZONTAL);
		maxIterBox = sfg::Box::Create(sfg::Box::Orientation::HORIZONTAL);

		xBox->Pack(xLabel); //Пихаем все в коробки
		xBox->Pack(xEntry);
		yBox->Pack(yLabel);
		yBox->Pack(yEntry);
		zoomBox->Pack(zoomLabel);
		zoomBox->Pack(zoomEntry);
		maxIterBox->Pack(maxIterLabel);
		maxIterBox->Pack(maxIterEntry);

		settingsBox->Pack(xBox); //пихаем коробки в основную коробку
		settingsBox->Pack(yBox);
		settingsBox->Pack(zoomBox);
		settingsBox->Pack(maxIterBox);
		settingsBox->Pack(previewButton);

		settingsBox->SetSpacing(5.f);
		settingsWindow->Add(settingsBox); //основную коробку в окно

										  //НАстройка окна с превью
		previewBox = sfg::Box::Create(sfg::Box::Orientation::VERTICAL);
		controlsLabel = sfg::Label::Create("Controls:\nWASD,\nShift,\n+/-");
		previewBox->Pack(controlsLabel);
		previewWindow->Add(previewBox);

		//Настройка окна экспорта
		widthLabel = sfg::Label::Create("width:"); //width
		widthEntry = sfg::Entry::Create("512");
		widthEntry->SetRequisition(sf::Vector2f(200.f, 10.f));
		heightLabel = sfg::Label::Create("height:"); //height
		heightEntry = sfg::Entry::Create("512");
		heightEntry->SetRequisition(sf::Vector2f(100.f, 10.f));
		SSAALabel = sfg::Label::Create("SSAA:");
		SSAAComboBox = sfg::ComboBox::Create();
		SSAAComboBox->AppendItem("None");
		SSAAComboBox->AppendItem("2x");
		SSAAComboBox->AppendItem("4x");
		SSAAComboBox->AppendItem("8x");
		SSAAComboBox->AppendItem("16x");
		SSAAComboBox->SelectItem(0);
		threadNumLabel = sfg::Label::Create("threads:");
		threadNumSpinbutton = sfg::SpinButton::Create(1.f, 64.f, 1.f);
		threadNumSpinbutton->SetRequisition(sf::Vector2f(80.f, 0.f));
		threadNumSpinbutton->SetValue(4.f);

		estimateTimeButton = sfg::Button::Create();
		estimateTimeButton->SetLabel("Estimate time");
		estimateTimeButton->GetSignal(sfg::Widget::OnLeftClick).Connect([this] {estimateTime(); });

		estimatedTimeLabel = sfg::Label::Create("Est. time: -");

		generateButton = sfg::Button::Create("Generate");
		generateButton->GetSignal(sfg::Widget::OnLeftClick).Connect([this] {generate(); });

		spinner = sfg::Spinner::Create();
		spinner->SetRequisition(sf::Vector2f(40.f, 40.f));

		generateTimeLabel = sfg::Label::Create("Not finished");

		colorLabel = sfg::Label::Create("Coloring:"); //Раскраска
		colorComboBox = sfg::ComboBox::Create();
		colorComboBox->AppendItem("Grayscale");
		colorComboBox->AppendItem("Hue cycle");

		filenameEntry = sfg::Entry::Create("filename");
		exportButton = sfg::Button::Create("Export");
		exportButton->GetSignal(sfg::Widget::OnLeftClick).Connect([this] {exportbut(); });

		autoMakeCheckButton = sfg::CheckButton::Create("Auto generate and export");
		zoomStepLabel = sfg::Label::Create("Zoom multipliyer");
		zoomStepEntry = sfg::Entry::Create();

		widthBox = sfg::Box::Create(sfg::Box::Orientation::HORIZONTAL);
		widthBox->Pack(widthLabel);
		widthBox->Pack(widthEntry);

		heigthBox = sfg::Box::Create(sfg::Box::Orientation::HORIZONTAL);
		heigthBox->Pack(heightLabel);
		heigthBox->Pack(heightEntry);

		threadBox = sfg::Box::Create(sfg::Box::Orientation::HORIZONTAL);
		threadBox->Pack(threadNumLabel);
		threadBox->Pack(threadNumSpinbutton);

		SSAABox = sfg::Box::Create(sfg::Box::Orientation::HORIZONTAL);
		SSAABox->Pack(SSAALabel);
		SSAABox->Pack(SSAAComboBox);

		exportBox = sfg::Box::Create(sfg::Box::Orientation::VERTICAL);
		exportBox->Pack(widthBox);
		exportBox->Pack(heigthBox);
		exportBox->Pack(SSAABox);
		exportBox->Pack(threadBox);
		exportBox->Pack(estimateTimeButton);
		exportBox->Pack(estimatedTimeLabel);
		exportBox->Pack(generateButton);
		exportBox->Pack(spinner);
		exportBox->Pack(generateTimeLabel);
		exportBox->Pack(colorLabel);
		exportBox->Pack(colorComboBox);
		exportBox->Pack(filenameEntry);
		exportBox->Pack(exportButton);
		exportBox->Pack(autoMakeCheckButton);
		exportBox->Pack(zoomStepLabel);
		exportBox->Pack(zoomStepEntry);

		exportBox->SetSpacing(5.f);

		exportWindow->Add(exportBox);

		//Добавляем окна на рабочий стол
		desktop.Add(settingsWindow);
		desktop.Add(previewWindow);
		desktop.Add(exportWindow);

		calcPreview();
	}
	void Loop()
	{
		sf::Clock clock;

		while (app_window.isOpen())
		{
			sf::Event event;
			while (app_window.pollEvent(event))
			{
				desktop.HandleEvent(event);

				if (event.type == sf::Event::Closed)
				{
					return;
				}

				if (event.type == sf::Event::MouseButtonPressed)
				{
					GetInputs();

					sf::Vector2i mousepos = sf::Mouse::getPosition(app_window);
					sf::Vector2f windpos = previewWindow->GetAbsolutePosition();



					int clickposx = sf::Mouse::getPosition(app_window).x - previewBox->GetChildren().at(1)->GetAbsolutePosition().x;
					int clickposy = sf::Mouse::getPosition(app_window).y - previewBox->GetChildren().at(1)->GetAbsolutePosition().y;
					if (clickposx > 0 && clickposx < 400 && clickposy > 0 && clickposy < 400)
					{
						inputOrbitX = (mpf_class)(clickposx - 200)/400 * inputZoom + inputX;
						inputOrbitY = -(mpf_class)(clickposy - 200)/400 * inputZoom + inputY;
					}
				}

				if (event.type == sf::Event::KeyPressed)
				{
					mpf_class step = 0.5;
					if (sf::Keyboard::isKeyPressed(sf::Keyboard::LShift)) // Управление
					{
						step = 0.15;
					}
					if (sf::Keyboard::isKeyPressed(sf::Keyboard::A))
					{
						inputX = inputX - step*inputZoom;
						SetInputs();
						calcPreview();
					}
					if (sf::Keyboard::isKeyPressed(sf::Keyboard::D))
					{
						inputX = inputX + step*inputZoom;
						SetInputs();
						calcPreview();
					}
					if (sf::Keyboard::isKeyPressed(sf::Keyboard::W))
					{
						inputY = inputY + step*inputZoom;
						SetInputs();
						calcPreview();
					}
					if (sf::Keyboard::isKeyPressed(sf::Keyboard::S))
					{
						inputY = inputY - step*inputZoom;
						SetInputs();
						calcPreview();
					}
					if (sf::Keyboard::isKeyPressed(sf::Keyboard::Add))
					{
						inputZoom = inputZoom * 0.5;
						SetInputs();
						calcPreview();
					}
					if (sf::Keyboard::isKeyPressed(sf::Keyboard::Subtract))
					{
						inputZoom = inputZoom * 2;
						SetInputs();
						calcPreview();
					}
				}

			}

			controlThreads();

			if (clock.getElapsedTime().asMicroseconds() >= 5000)
			{
				// Update() takes the elapsed time in seconds.
				desktop.Update(static_cast<float>(clock.getElapsedTime().asMicroseconds()) / 1000000.f);
				clock.restart();
			}
			app_window.clear();
			sfgui.Display(app_window);
			app_window.display();

			//std::cout << inputOrbitX << ", "  << inputOrbitY << std::endl;
			//std::cout << previewBox->GetChildren().at(0)->GetAbsolutePosition().x << ", " << previewBox->GetChildren().at(0)->GetAbsolutePosition().y << std::endl;
			/*std::cout << "export: " << exportWindow->GetRequisition().x << ", " << exportWindow->GetRequisition().y << "\n";
			std::cout << "settings: " << settingsWindow->GetRequisition().x << ", " << settingsWindow->GetRequisition().y << "\n";
			std::cout << "preview: " << previewWindow->GetRequisition().x << ", " << previewWindow->GetRequisition().y << "\n";*/
		}
	}
};

int main() {


	GUI app;
	app.SetupGUI();
	app.Loop();

	/*sf::RenderWindow window(sf::VideoMode(800, 600), "fdf");
	window.setFramerateLimit(30);

	Set mandelbrot(mpf_class(-1), mpf_class(0), 4, 360, 1, 800, 600);
	mandelbrot.Generate(4);

	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
			{
				return 0;
			}

		}

		window.clear();

		sf::Image pic = mandelbrot.GetImage(HUECYCLE);

		sf::Texture text;
		text.loadFromImage(pic);
		sf::Sprite s;
		s.setTexture(text);

		window.draw(s);
		window.display();
	}*/


	return EXIT_SUCCESS;
}