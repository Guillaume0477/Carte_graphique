
#include <cassert>
#include <cstdio>
#include <vector>
#include <set>
#include <string>

#define GLEW_NO_GLU
#include "GL/glew.h"

#include "window.h"

#ifndef _MSC_VER
  #define GK_CALLBACK
#else
  #define GK_CALLBACK __stdcall
#endif


static int width;
static int height;
int window_width( )
{
    return width;
}
int window_height( )
{
    return height;
}

static std::vector<unsigned char> key_states;
int key_state( const SDL_Keycode key )
{
    SDL_Scancode code= SDL_GetScancodeFromKey(key);
    assert((size_t) code < key_states.size());
    return (int)  key_states[code];
}
void clear_key_state( const SDL_Keycode key )
{
    SDL_Scancode code= SDL_GetScancodeFromKey(key);
    assert((size_t) code < key_states.size());
    key_states[code]= 0;
}

static SDL_KeyboardEvent last_key;
SDL_KeyboardEvent key_event( )
{
    return last_key;
}
void clear_key_event( )
{
    last_key.type= 0;
    last_key.keysym.sym= 0;
}

static SDL_TextInputEvent last_text;
SDL_TextInputEvent text_event( )
{
    return last_text;
}
void clear_text_event( )
{
    last_text.text[0]= 0;
}

static std::string last_drop;
const char *drop_event( )
{
    return last_drop.c_str();
}
void clear_drop_event( )
{
    last_drop.clear();
}

static SDL_MouseButtonEvent last_button;
SDL_MouseButtonEvent button_event( )
{
    return last_button;
}
void clear_button_event( )
{
    last_button.state= 0;
}

static SDL_MouseWheelEvent last_wheel;
SDL_MouseWheelEvent wheel_event( )
{
    return last_wheel;
}
void clear_wheel_event( )
{
    last_wheel.x= 0;
    last_wheel.y= 0;
}


//! boucle de gestion des evenements de l'application.
int draw( );    //!< declaration anticipee.

int run( Window window )
{
    SDL_StartTextInput();

    SDL_EventState(SDL_DROPFILE, SDL_ENABLE);
    last_drop= std::string();
    
    int stop= 0;
    while(stop == 0)
    {
        // gestion des evenements
        SDL_Event event;
        while(SDL_PollEvent(&event))
        {
            switch(event.type)
            {
                case SDL_WINDOWEVENT:
                    // redimensionner la fenetre...
                    if(event.window.event == SDL_WINDOWEVENT_RESIZED)
                    {
                        // conserve les dimensions de la fenetre
                        width= event.window.data1;
                        height= event.window.data2;
                        SDL_SetWindowSize(window, width, height);

                        // ... et le viewport opengl
                        glViewport(0, 0, width, height);
                    }
                    break;

                case SDL_DROPFILE:
                    last_drop.assign(event.drop.file);
                    SDL_free(event.drop.file);
                    break;
                
                case SDL_TEXTINPUT:
                    // conserver le dernier caractere
                    last_text= event.text;
                    break;
                
                case SDL_KEYDOWN:
                    // modifier l'etat du clavier
                    if((size_t) event.key.keysym.scancode < key_states.size())
                    {
                        key_states[event.key.keysym.scancode]= 1;
                        last_key= event.key;    // conserver le dernier evenement
                    }
                    
                    // fermer l'application
                    if(event.key.keysym.sym == SDLK_ESCAPE)
                        stop= 1;
                    break;
                
                case SDL_KEYUP:
                    // modifier l'etat du clavier
                    if((size_t) event.key.keysym.scancode < key_states.size())
                    {
                        key_states[event.key.keysym.scancode]= 0;
                        last_key= event.key;    // conserver le dernier evenement
                    }
                    break;

                case SDL_MOUSEBUTTONDOWN:
                case SDL_MOUSEBUTTONUP:
                    last_button= event.button;
                    break;
                
                case SDL_MOUSEWHEEL:
                    last_wheel= event.wheel;
                    break;
                
                case SDL_QUIT:
                    stop= 1;            // fermer l'application
                    break;
            }
        }

        // dessiner
        if(draw() < 1)
            stop= 0;    // fermer l'application si draw() renvoie 0 ou -1...
        
        // presenter le resultat
        SDL_GL_SwapWindow(window);
    }

    SDL_StopTextInput();
    return 0;
}

//! creation d'une fenetre pour l'application.
Window create_window( const int w, const int h )
{
    // init sdl
    if(SDL_Init(SDL_INIT_EVERYTHING | SDL_INIT_NOPARACHUTE ) < 0)
    {
        printf("[error] SDL_Init() failed:\n%s\n", SDL_GetError());
        return NULL;
    }

    // enregistre le destructeur de sdl
    atexit(SDL_Quit);

    // creer la fenetre 
    Window window= SDL_CreateWindow("gKit",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, w, h,
        SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
    if(window == NULL)
    {
        printf("[error] SDL_CreateWindow() failed.\n");
        return NULL;
    }

    // recupere l'etat du clavier
    int keys;
    const unsigned char *state= SDL_GetKeyboardState(&keys);
    key_states.assign(state, state + keys);
    
    // conserve les dimensions de la fenetre
    width= w;
    height= h;
    
    SDL_SetWindowDisplayMode(window, NULL);
    return window;
}

void release_window( Window window )
{
    SDL_DestroyWindow(window);
}


//! affiche les messages d'erreur opengl. (contexte debug core profile necessaire).
static
void GK_CALLBACK debug( GLenum source, GLenum type, unsigned int id, GLenum severity, GLsizei length,
    const char *message, void *userParam )
{
    static std::set<std::string> log;
    if(log.insert(message).second == false)
        // le message a deja ete affiche, pas la peine de recommencer 60 fois par seconde.
        return;

    if(severity == GL_DEBUG_SEVERITY_HIGH)
        printf("[openGL error]\n%s\n", message);
#ifndef GK_RELEASE
    else if(severity == GL_DEBUG_SEVERITY_MEDIUM)
        printf("[openGL warning]\n%s\n", message);
    else
        printf("[openGL message]\n%s\n", message);
#endif
}

//! cree et configure un contexte opengl
Context create_context( Window window, const int major, const int minor )
{
    if(window == NULL)
        return NULL;

    // configure la creation du contexte opengl core profile, debug profile
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, major);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, minor);
#ifndef GK_RELEASE    
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_DEBUG_FLAG);
#endif
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 15);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

    Context context= SDL_GL_CreateContext(window);
    if(context == NULL)
    {
        printf("[error] creating openGL context.\n");
        return NULL;
    }

    SDL_GL_SetSwapInterval(1);

    // initialise les extensions opengl
    glewExperimental= 1;
    GLenum err= glewInit();
    if(err != GLEW_OK)
    {
        printf("[error] loading extensions\n%s\n", glewGetErrorString(err));
        SDL_GL_DeleteContext(context);
        return NULL;
    }

    // purge les erreurs opengl generees par glew !
    while(glGetError() != GL_NO_ERROR) {;}

#ifndef GK_RELEASE        
    // configure l'affichage des messages d'erreurs opengl, si l'extension est disponible
    if(GLEW_ARB_debug_output)
    {
        glDebugMessageControlARB(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, NULL, GL_TRUE);
        glDebugMessageCallbackARB(debug, NULL);
        glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS_ARB);
    }
#endif
    
    return context;
}

void release_context( Context context )
{
    SDL_GL_DeleteContext(context);
}

